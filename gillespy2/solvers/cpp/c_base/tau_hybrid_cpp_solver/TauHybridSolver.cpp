/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2021 GillesPy2 developers.
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

namespace Gillespy
{
	namespace TauHybrid
	{
		template <typename T>
		using DelayedExecutionQueue = std::priority_queue<EventExecution, std::vector<EventExecution>, T>;

		bool interrupted = false;

		void signalHandler(int signum)
		{
			interrupted = true;
		}

		void TauHybridCSolver(HybridSimulation *simulation, std::vector<Event> &events, const double tau_tol)
		{
			if (simulation == NULL)
			{
				return;
			}

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
			} else
			{
				y = y0;
			}
			Integrator sol(simulation, y, GPY_HYBRID_RELTOL, GPY_HYBRID_ABSTOL);

			// Tau selector initialization. Used to select a valid tau step.
			TauArgs<double> tau_args = initialize(model, tau_tol);

			// Simulate for each trajectory
			for (int traj = 0; traj < num_trajectories; traj++)
			{
				if (traj > 0)
				{
					sol.reinitialize(y0);
				}

				// Initialize each species with their respective user modes.
				for (int spec_i = 0; spec_i < num_species; ++spec_i)
				{
					HybridSpecies *spec = &simulation->species_state[spec_i];
					spec->partition_mode = spec->user_mode == SimulationState::DYNAMIC
										   ? SimulationState::DISCRETE
										   : spec->user_mode;
					simulation->trajectories[traj][0][spec_i] = spec->base_species->initial_population;
				}

				// Population/concentration state values for each species.
				// TODO: change back double -> hybrid_state, once we figure out how that works
				std::vector<double> current_state(num_species);
				std::vector<int> current_populations(num_species);

				// Initialize the species population for the trajectory.
				for (int spec_i = 0; spec_i < num_species; ++spec_i)
				{
					current_state[spec_i] = species[spec_i].initial_population;
					current_populations[spec_i] = species[spec_i].initial_population;
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
				// This is a temporary fix. Ideally, invalid state should allow for integrator options change.
				// For now, a "guard" is put in place to prevent potentially infinite loops from occurring.
				unsigned int integration_guard = 1000;

				// delay_queue:
				// Delayed execution objects which are guaranteed to (eventually) fire at its time.
				// Once placed in this queue, an execution is only removed when it's ready to fire.
				DelayedExecutionQueue<std::greater<EventExecution>> delay_queue;

				// volatile_queue:
				// Delayed execution objects which must be re-evaluated at each root.
				// If any volatile executions transition from True -> False, they are removed.
				// Otherwise, they are treated just like a delay_queue, by triggering when ready.
				std::list<EventExecution> volatile_queue;

				while (integration_guard > 0 && simulation->current_time < simulation->end_time)
				{
					// Compute current propensity values based on existing state.
					for (unsigned int rxn_i = 0; rxn_i < num_reactions; ++rxn_i)
					{
						HybridReaction &rxn = simulation->reaction_state[rxn_i];
						double propensity = 0.0;
						switch (rxn.mode)
						{
						case SimulationState::CONTINUOUS:
							propensity = Reaction::propensity(rxn_i, current_state.data());
							break;
						case SimulationState::DISCRETE:
							propensity = Reaction::propensity(rxn_i, current_populations.data());
							break;
						default:
							break;
						}
						sol.data.propensities[rxn_i] = propensity;
					}

					// Expected tau step is determined.
					tau_step = select(
							model,
							tau_args,
							tau_tol,
							simulation->current_time,
							save_time,
							sol.data.propensities,
							current_populations
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
					next_time = events.empty()
						? simulation->current_time + tau_step
						: simulation->current_time + std::min(tau_step, increment);

					// The integration loop continues until a valid solution is found.
					// Any invalid Tau steps (which cause negative populations) are discarded.
					sol.save_state();
					do
					{
						// Integration Step
						// For deterministic reactions, the concentrations are updated directly.
						// For stochastic reactions, integration updates the rxn_offsets vector.
						IntegrationResults result = sol.integrate(&next_time);
						if (sol.status == IntegrationStatus::BAD_STEP_SIZE)
						{
							invalid_state = true;
							// Breaking early causes `invalid_state` to remain set,
							//   resulting in an early termination of the trajectory.
							break;
						}

						// The integrator has, at this point, been validated.
						// Any errors beyond this point is assumed to be a stochastic state failure.
						invalid_state = false;

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
								while (rxn_state >= 0)
								{
									// "Fire" a reaction by recording changes in dependent species.
									// If a negative value is detected, break without saving changes.
									for (int spec_i = 0; spec_i < num_species; ++spec_i)
									{
										population_changes[spec_i] +=
												model.reactions[rxn_i].species_change[spec_i];
										if (current_state[spec_i] + population_changes[spec_i] < 0)
										{
											invalid_state = true;
										}
									}

									rxn_state += log(urn.next());
								}
								result.reactions[rxn_i] = rxn_state;
								break;

							case SimulationState::CONTINUOUS:
							default:
								break;
							}
						}

						// Positive reaction state means a negative population was detected.
						// Only update state with the given population changes if valid.
						if (invalid_state)
						{
							sol.restore_state();
							tau_step /= 2;
							next_time = simulation->current_time + tau_step;
						} else
						{
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
								current_populations[p_i] = (int) current_state[p_i];
							}
						}
					} while (invalid_state);

					// Invalid state after the do-while loop implies that an unrecoverable error has occurred.
					// While prior results are considered usable, the current integration results are not.
					// Calling `continue` with an invalid state will discard the results and terminate the trajectory.
					integration_guard = invalid_state
										? integration_guard - 1
										: 1000;

					if (!events.empty())
					{
						double t = simulation->current_time;
						double *event_state = N_VGetArrayPointer(sol.y);
						auto compare = [t, event_state](EventExecution &lhs, EventExecution &rhs) -> bool
						{
							return lhs.priority(t, event_state) < rhs.priority(t, event_state);
						};
						std::priority_queue<EventExecution, std::vector<EventExecution>, decltype(compare)>
								trigger_queue(compare);

						for (auto &event : events)
						{
							if (event.trigger(t, event_state))
							{
								double delay = event.delay(t, event_state);
								if (delay <= 0)
								{
									// Immediately put EventExecution on "triggered" pile
									trigger_queue.emplace(event.get_execution(t, event_state, num_species));
								}
								else if (event.is_persistent())
								{
									// Put EventExecution on "delayed" pile
									delay_queue.emplace(event.get_execution(delay, event_state, num_species));
								}
								else
								{
									// Delayed, but must be re-checked on every iteration.
									volatile_queue.emplace_back(event.get_execution(delay, event_state, num_species));
								}
							}
						}

						for (auto vol_event = volatile_queue.begin(); vol_event != volatile_queue.end(); ++vol_event)
						{
							// Execution objects in the volatile queue must remain True until execution.
							// Remove any execution objects which transitioned to False before execution.
							if (!vol_event->trigger(t, event_state))
							{
								vol_event = volatile_queue.erase(vol_event);
							}
							else if (vol_event->get_execution_time() >= t)
							{
								trigger_queue.push(*vol_event);
								vol_event = volatile_queue.erase(vol_event);
							}
						}

						while (!delay_queue.empty())
						{
							auto &event = delay_queue.top();
							if (event.get_execution_time() < t)
							{
								break;
							}
							trigger_queue.push(event);
							delay_queue.pop();
						}

						while (!trigger_queue.empty())
						{
							auto event = trigger_queue.top();
							event.execute(t, event_state);
							trigger_queue.pop();
						}
					}

					// Output the results for this time step.
					sol.refresh_state();
					simulation->current_time = next_time;

					// Seek forward, writing out any values on the timeline which are on current timestep range.
					while (save_idx < simulation->number_timesteps && save_time <= next_time)
					{
						for (int spec_i = 0; spec_i < num_species; ++spec_i)
						{
							simulation->trajectories[traj][save_idx][spec_i] = current_state[spec_i];
						}
						save_time = simulation->timeline[++save_idx];
					}
				}

				if (integration_guard == 0)
				{
					std::cerr
							<< "[Trajectory #" << traj << "] "
							<< "Integration guard triggered; problem space too stiff at t="
							<< simulation->current_time << std::endl;
				}

				// End of trajectory
				delete[] population_changes;
			}
		}
	}
}
