/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2022 GillesPy2 developers.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <csignal> //Included for timeout signal handling
#include <random>
#include <queue>
#include <list>
#include "cvode.h" // prototypes for CVODE fcts., consts.
#include "nvector_serial.h"  // access to serial N_Vector
#include "sunlinsol_spgmr.h"  //access to SPGMR SUNLinearSolver
#include "cvode_spils.h" // access to CVSpils interface
#include "sundials_types.h"  // defs. of realtype, sunindextype
#include "sundials_math.h"  // contains the macros ABS, SUNSQR, EXP
#include "TauHybridSolver.h"
#include "HybridModel.h"
#include "integrator.h"
#include "tau.h"

static void silent_error_handler(int error_code, const char *module, const char *function_name,
						  char *message, void *eh_data);

#define INTEGRATION_GUARD_MAX 1

namespace Gillespy
{
	static volatile bool interrupted = false;

	GPY_INTERRUPT_HANDLER(signal_handler, {
		interrupted = true;
	})

	namespace TauHybrid
	{
        bool TakeIntegrationStep(Integrator&sol, double next_time, int*population_changes, std::vector<double>*current_state, std::set<unsigned int>*rxn_roots,  HybridSimulation*simulation, URNGenerator*urn){
                TauHybrid::TakeIntegrationStep(sol, nextime, population_changes, current_state, rxn_roots, simulation, urn, -1);
        }

        bool TakeIntegrationStep(Integrator&sol, double next_time, int*population_changes, std::vector<double>*current_state, std::set<unsigned int>*rxn_roots,  HybridSimulation*simulation, URNGenerator*urn, int only_reaction_to_fire){
            Model<double> &model = *(simulation->model);
            int num_species = model.number_species;
            // Integration Step
            // For deterministic reactions, the concentrations are updated directly.
            // For stochastic reactions, integration updates the rxn_offsets vector.
            IntegrationResults result = sol.integrate(&next_time, event_roots, rxn_roots, urn);
            if (sol.status == IntegrationStatus::BAD_STEP_SIZE)
            {
                std::cerr << "IntegrationStatus::BAD_STEP_SIZE sol.status="<<sol.status<<"\n";
                return false;
            }
            else
            {
                // The integrator has, at this point, been validated.
                // Any errors beyond this point is assumed to be a stochastic state failure.

                // 0-initialize our population_changes array.
                for (int p_i = 0; p_i < num_species; ++p_i)
                {
                    population_changes[p_i] = 0;
                }

                // Start with the species concentration as a baseline value.
                // Stochastic reactions will update populations relative to their concentrations.
                for (int spec_i = 0; spec_i < num_species; ++spec_i)
                {
                    current_state[spec_i] = result.concentrations[spec_i]; 
                }

                if (!rxn_roots.empty())
                {
                    // "Direct" roots found; these are executed manually
                    for (unsigned int rxn_i : rxn_roots)
                    {
                        // "Fire" a reaction by recording changes in dependent species.
                        // If a negative value is detected, break without saving changes.
                        for (int spec_i = 0; spec_i < num_species; ++spec_i)
                        {
                            // Unlike the Tau-leaping version of reaction firings,
                            // it is not possible to have a negative state occur in direct reactions.
                            population_changes[spec_i] += model.reactions[rxn_i].species_change[spec_i];
                            result.reactions[rxn_i] = log(urn.next());
                        }
                    }
                    rxn_roots.clear();
                }
                else
                {
                    // The newly-updated reaction_states vector may need to be reconciled now.
                    // A positive reaction_state means reactions have potentially fired.
                    // NOTE: it is possible for a population to swing negative, where a smaller Tau is needed.
                    for (int rxn_i = 0; rxn_i < num_reactions; ++rxn_i)
                    {
                        // Temporary variable for the reaction's state.
                        // Does not get updated unless the changes are deemed valid.
                        double rxn_state = result.reactions[rxn_i];

                        switch (simulation->reaction_state[rxn_i].mode)
                        {
                        case SimulationState::DISCRETE:
                            unsigned int rxn_count = 0;
                            if(only_reaction_to_fire == rxn_i){
                                    std::cerr << "Firing single SSA reaction "<< rxn_i<<"\n";
                                    rxn_state = log(urn.next());
                                
                            }else if(rxn_state > 0){
                                while (rxn_state >= 0) {
                                    // "Fire" a reaction by recording changes in dependent species.
                                    // If a negative value is detected, break without saving changes. (<-- this is not implemented)
                                    rxn_state += log(urn.next());
                                    rxn_count++;
                                }
                            }
                            if(rxn_count > 0){
                                for (int spec_i = 0; spec_i < num_species; ++spec_i) {
                                    population_changes[spec_i] += model.reactions[rxn_i].species_change[spec_i] * rxn_count;
                                }
                            }
                            result.reactions[rxn_i] = rxn_state;
                            break;

                        case SimulationState::CONTINUOUS:
                        default:
                            break;
                        }
                    }
                }
        }

        bool NegativeStateCheck(int num_species, int*population_changes, std::vector<double>*current_state){
            // Explicitly check for invalid population state, now that changes have been tallied.
            bool invalid_state = false;
            for (int spec_i = 0; spec_i < num_species; ++spec_i)
            {
                if (current_state[spec_i] + population_changes[spec_i] < 0)
                {
                    invalid_state = true;
                    break;
                }
            }
            return invalid_state;
        }

		void TauHybridCSolver(
				HybridSimulation *simulation,
				std::vector<Event> &events,
				Logger &logger,
				double tau_tol,
				SolverConfiguration config)
		{
			if (simulation == NULL)
			{
				return;
			}
			GPY_INTERRUPT_INSTALL_HANDLER(signal_handler);

			Model<double> &model = *(simulation->model);
			int num_species = model.number_species;
			int num_reactions = model.number_reactions;
			int num_trajectories = simulation->number_trajectories;
			std::unique_ptr<Species<double>[]> &species = model.species;
			double increment = simulation->timeline[1] - simulation->timeline[0];

			URNGenerator urn(simulation->random_seed);
			// The contents of y0 are "stolen" by the integrator.
			// Do not attempt to directly use y0 after being passed to sol!
			N_Vector y0 = init_model_vector(model, urn);
			N_Vector y;
			if (num_trajectories > 0)
			{
				y = init_model_vector(model, urn);
			}
			else
			{
				y = y0;
			}
			Integrator sol(simulation, y, config.rel_tol, config.abs_tol);
			if (logger.get_log_level() == LogLevel::CRIT)
			{
				sol.set_error_handler(silent_error_handler);
			}

			// Configure user-specified solver tolerances.
			if (!sol.configure(config))
			{
				logger.warn() << "Received invalid tolerances: {"
					<< "rtol = " << config.rel_tol
					<< ", atol = " << config.abs_tol
					<< ", max_step = " << config.max_step
					<< "}" << std::endl;
			}

			// Tau selector initialization. Used to select a valid tau step.
			TauArgs<double> tau_args = initialize(model, tau_tol);

			// Simulate for each trajectory
			for (int traj = 0; !interrupted && traj < num_trajectories; traj++)
			{
				if (traj > 0)
				{
					sol.reinitialize(y0);
				}

				// Population/concentration state values for each species.
				// TODO: change back double -> hybrid_state, once we figure out how that works
				EventList event_list;
				std::vector<double> current_state(num_species);

				// Initialize the species population for the trajectory.
				for (int spec_i = 0; spec_i < num_species; ++spec_i)
				{
					current_state[spec_i] = species[spec_i].initial_population;
					simulation->current_state[spec_i] = current_state[spec_i];
				}
				simulation->reset_output_buffer(traj);
				simulation->output_buffer_range(std::cout);

				// Check for initial event triggers at t=0 (based on initial_value of trigger)
				std::set<int> event_roots;
				std::set<unsigned int> rxn_roots;
				if (event_list.evaluate_triggers(current_state.data(), simulation->current_time))
				{
					double *event_state = N_VGetArrayPointer(sol.y);
					event_list.evaluate(current_state.data(), num_species, simulation->current_time, event_roots);
					std::copy(current_state.begin(), current_state.end(), event_state);
					sol.refresh_state();
				}

				// Initialize each species with their respective user modes.
				for (int spec_i = 0; spec_i < num_species; ++spec_i)
				{
					HybridSpecies *spec = &simulation->species_state[spec_i];
					spec->partition_mode = spec->user_mode == SimulationState::DYNAMIC
										   ? SimulationState::DISCRETE
										   : spec->user_mode;
				}

				// SIMULATION STEP LOOP
				int save_idx = 1;
				double next_time;
				double tau_step = 0.0;
				double save_time = simulation->timeline[save_idx];

				// Temporary array to store changes to dependent species.
				// Should be 0-initialized each time it's used.
				int *population_changes = new int[num_species];
				simulation->current_time = 0;

				// An invalid simulation state indicates that an unrecoverable error has occurred,
				//   and the trajectory should terminate early.
				bool invalid_state = false;


				while (!interrupted && !invalid_state && simulation->current_time < simulation->end_time)
				{
					// Compute current propensity values based on existing state.
					for (unsigned int rxn_i = 0; rxn_i < num_reactions; ++rxn_i)
					{
						HybridReaction &rxn = simulation->reaction_state[rxn_i];
						sol.data.propensities[rxn_i] = rxn.ssa_propensity(current_state.data());
					}
					if (interrupted)
						break;

					// Expected tau step is determined.
					tau_step = select<double, double>(
							model,
							tau_args,
							tau_tol,
							simulation->current_time,
							save_time,
							sol.data.propensities,
							current_state
					);
					partition_species(
							simulation->reaction_state,
							simulation->species_state,
							sol.data.propensities,
							current_state,
							tau_step,
							tau_args
					);
					flag_det_rxns(
							simulation->reaction_state,
							simulation->species_state
					);
					update_species_state(simulation->species_state, current_state);
					create_differential_equations(simulation->species_state, simulation->reaction_state);

					// Determine what the next time point is.
					// This will become current_time on the next iteration.
					// If a retry with a smaller tau_step is deemed necessary, this will change.
					next_time = simulation->current_time + tau_step;

					// Ensure that any previous changes to the current state is reflected by the integrator.
					std::copy(current_state.begin(), current_state.end(), N_VGetArrayPointer(sol.y));

					// The integration loop continues until a valid solution is found.
					// Any invalid Tau steps (which cause negative populations) are discarded.
					sol.save_state();

					// This is a temporary fix. Ideally, invalid state should allow for integrator options change.
					// For now, a "guard" is put in place to prevent potentially infinite loops from occurring.
					unsigned int integration_guard = INTEGRATION_GUARD_MAX;

					do
					{
                        invalid_state = TauHybrid::TakeIntegrationStep(sol, next_time, populaton_changes, current_state, rxn_roots, simulation, urn);
                        if(!invalid_state){
                            invalid_state = TauHybrid::NegativeStateCheck(num_species, current_state, population_changes);
                        }

                        // If state is invalid, we took too agressive tau step and need to take a single SSA step forward
						if (invalid_state) {
                            // Re-Initialize the species population for the trajectory.
                            for (int spec_i = 0; spec_i < num_species; ++spec_i) {
                                current_state[spec_i] = species[spec_i].initial_population;
                                simulation->current_state[spec_i] = current_state[spec_i];
                            }
                            // Restore the solver to the intial step state
                            sol.restore_state();
                            //The above code blocks should re-initialize the state, is this correct?

                            // estimate the time to the first stochatic reaction by assuming constant propensities
                            double min_tau = 0.0;
                            unsigned int rxn_selected = -1
                            for (unsigned int rxn_i = 0; rxn_i < num_reactions; ++rxn_i) {
                                HybridReaction &rxn = simulation->reaction_state[rxn_i];
                                //estimate the zero crossing time

                                // Python code:
                                //rxn_times[rname] = -1* curr_state[rname] / propensities[rname]
                                // C++ code:
                                double *rxn_value = sol.get_reaction_state()
                                double est_tau = -1*  rxn_value[rxn_i] / rxn.ssa_propensity(current_state.data());

                                if(rxn_selected == -1 || est_tau < min_tau ){
                                    min_tau = est_tau;
                                    rxn_selected = rxn_i;
                                }
                            }
                            if(rxn_selected == -1){
                                invalid_state = true;
                                break;
                            }
                            // Use the found tau-step for single SSA
					        next_time = simulation->current_time + min_tau;

                            // Integreate the system forward
                            invalid_state = TauHybrid::TakeIntegrationStep(sol, next_time, populaton_changes, current_state, rxn_roots, simulation, urn, rxn_selected);
                            if(!invalid_state){
                                invalid_state = TauHybrid::NegativeStateCheck(num_species, current_state, population_changes);
                            }
                        }

						// Positive reaction state means a negative population was detected.
						// Only update state with the given population changes if valid.
						if (invalid_state) {
                            //Got an invalid state after the SSA step
                            break;

						} else {
                            // "Permanently" update the rxn_state and populations.
							for (int p_i = 0; p_i < num_species; ++p_i)
							{
								if (!simulation->species_state[p_i].boundary_condition)
								{
									// Boundary conditions are not modified directly by reactions.
									// As such, population dx in stochastic regime is not considered.
									// For deterministic species, their effective dy/dt should always be 0.
									current_state[p_i] += population_changes[p_i];
									result.concentrations[p_i] = current_state[p_i];
								}
							}
						}
					} while (invalid_state && !interrupted);

					if (interrupted)
						break;
					else if (invalid_state)
					{
						// Invalid state after the do-while loop implies that an unrecoverable error has occurred.
						// While prior results are considered usable, the current integration results are not.
						// Calling `continue` with an invalid state will discard the results and terminate the trajectory.
						logger.err()
								<< "[Trajectory #" << traj << "] "
								<< "Integration guard triggered; problem space too stiff at t="
								<< simulation->current_time << std::endl;
						simulation->set_status(HybridSimulation::LOOP_OVER_INTEGRATE);
						continue;
					}
					else
					{
						integration_guard = INTEGRATION_GUARD_MAX;
					}

					// ===== <EVENT HANDLING> =====
					if (!event_list.has_active_events())
					{
						if (event_list.evaluate_triggers(N_VGetArrayPointer(sol.y), next_time))
						{
							sol.restore_state();
							sol.use_events(events, simulation->reaction_state);
							sol.enable_root_finder();
							continue;
						}
					}
					else
					{
						double *event_state = N_VGetArrayPointer(sol.y);
						if (!event_list.evaluate(event_state, num_species, next_time, event_roots))
						{
							sol.disable_root_finder();
						}
						std::copy(event_state, event_state + num_species, current_state.begin());
					}
					// ===== </EVENT HANDLING> =====

					// Output the results for this time step.
					sol.refresh_state();
					simulation->current_time = next_time;

					// Seek forward, writing out any values on the timeline which are on current timestep range.
					while (save_idx < simulation->number_timesteps && save_time <= next_time)
					{
						for (int spec_i = 0; spec_i < num_species; ++spec_i)
						{
							simulation->current_state[spec_i] = current_state[spec_i];
						}
						simulation->output_buffer_range(std::cout, save_idx++);
						save_time = simulation->timeline[save_idx];
					}
				}

				// End of trajectory
				delete[] population_changes;
			}

			if (interrupted)
			{
				simulation->set_status(HybridSimulation::OK);
			}
		}
	}
}

void silent_error_handler(int error_code, const char *module, const char *function_name,
						  char *message, void *eh_data)
{
	// Do nothing
}
