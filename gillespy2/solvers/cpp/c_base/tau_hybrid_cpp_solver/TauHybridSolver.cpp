#include <iostream>
#include <csignal> //Included for timeout signal handling
#include <random>
#include "cvode.h" // prototypes for CVODE fcts., consts.
#include "nvector_serial.h"  // access to serial N_Vector
#include "sunlinsol_spgmr.h"  //access to SPGMR SUNLinearSolver
#include "cvode_spils.h" // access to CVSpils interface
#include "sundials_types.h"  // defs. of realtype, sunindextype
#include "sundials_math.h"  // contains the macros ABS, SUNSQR, EXP
#include "TauHybridSolver.h"
#include "HybridModel.h"
#include "integrator.h"
// #include "statistics.h"
#include "tau.h"
using namespace Gillespy;

namespace Gillespy::TauHybrid {
	bool interrupted = false;

	void signalHandler(int signum)
	{
		interrupted = true;
	}

	void TauHybridCSolver(HybridSimulation *simulation, const double tau_tol)
	{
		if (simulation == NULL) {
			return;
		}

		Model &model = *(simulation->model);
		int num_species = model.number_species;
		int num_reactions = model.number_reactions;
		int num_trajectories = simulation->number_trajectories;
		std::unique_ptr<Species[]> &species = model.species;
		double increment = simulation->timeline[1] - simulation->timeline[0];

		// Tau selector initialization. Used to select a valid tau step.
		TauArgs tau_args = initialize(model, tau_tol);

		//copy initial state for each trajectory
		for(int s = 0; s < num_species; s++){
			simulation->trajectoriesODE[0][0][s] = species[s].initial_population;
		}

		//Simulate for each trajectory
		for(int traj = 0; traj < num_trajectories; traj++){
			// Population/concentration state values for each species.
			// TODO: change back double -> hybrid_state, once we figure out how that works
			std::vector<double> current_state(num_species);
			std::vector<int> current_populations(num_species);

			URNGenerator urn;
			// The contents of y0 are "stolen" by the integrator.
			// Do not attempt to directly use y0 after being passed to sol!
			N_Vector y0 = init_model_vector(model, urn);
			Integrator sol(simulation, y0, GPY_HYBRID_RELTOL, GPY_HYBRID_ABSTOL);

			// Initialize the species population for the trajectory.
			for (int spec_i = 0; spec_i < num_species; ++spec_i) {
				current_state[spec_i] = species[spec_i].initial_population;
				current_populations[spec_i] = species[spec_i].initial_population;
			}

			// Represents the largest valid index of the output vector(s), y (and y0).
			int rxn_offset_boundary  = num_species + num_reactions;

			// SIMULATION STEP LOOP
			double next_time;
			double tau_step = 0.0;
			double save_time = simulation->timeline[1];

			// Temporary array to store changes to dependent species.
			// Should be 0-initialized each time it's used.
			int *population_changes = new int[num_species];
			simulation->current_time = 0;
			while (simulation->current_time < simulation->end_time) {
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

				// 0-initialize our population_changes array.
				for (int p_i = 0; p_i < num_species; ++p_i) {
					population_changes[p_i] = 0;
				}

				// Determine what the next time point is.
				// This will become current_time on the next iteration.
				// If a retry with a smaller tau_step is deemed necessary, this will change.
				next_time = simulation->current_time + tau_step;

				// Integration Step
				// For deterministic reactions, the concentrations are updated directly.
				// For stochastic reactions, integration updates the rxn_offsets vector.
				// flag = CVode(cvode_mem, next_time, y0, &next_time, CV_NORMAL);
				IntegrationResults result = sol.integrate(&next_time);
				bool invalid_state = false;
				do {
					// Start with the species concentration as a baseline value.
					// Stochastic reactions will update populations relative to their concentrations.
					for (int spec_i = 0; spec_i < num_species; ++spec_i) {
						current_state[spec_i] = result.concentrations[spec_i];
					}

					// The newly-updated reaction_states vector may need to be reconciled now.
					// A positive reaction_state means reactions have potentially fired.
					// NOTE: it is possible for a population to swing negative, where a smaller Tau is needed.
					for (int rxn_i = 0; rxn_i < num_reactions; ++rxn_i) {
						// Temporary variable for the reaction's state.
						// Does not get updated unless the changes are deemed valid.
						double rxn_state = result.reactions[rxn_i];

						switch (sol.data.reaction_state[rxn_i].mode) {
						case SimulationState::DISCRETE:
							while (rxn_state >= 0) {
								// "Fire" a reaction by recording changes in dependent species.
								// If a negative value is detected, break without saving changes.
								for (int spec_i = 0; spec_i < num_species; ++spec_i) {
									population_changes[spec_i] +=
										model.reactions[rxn_i].species_change[spec_i];
									if (current_state[spec_i] + population_changes[spec_i] < 0) {
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
					if (invalid_state) {
						// TODO: invalidate reaction timestep, reset integrator state, try again
						std::cerr << "Integration step failed; RIP" << std::endl;
						invalid_state = false;
					}
					else {
						// "Permanently" update the rxn_state and populations.
						for (int p_i = 0; p_i < num_species; ++p_i) {
							current_state[p_i] += population_changes[p_i];
							current_populations[p_i] = current_state[p_i];
							result.concentrations[p_i] = current_state[p_i];
						}
					}
				} while (invalid_state);

				// Output the results for this time step.
				sol.reset(next_time);
				simulation->current_time = next_time;
				
				while (save_time <= next_time) {
					// Write each species, one at a time (from ODE solution)
					for (int spec_i = 0; spec_i < num_species; ++spec_i) {
						simulation->trajectoriesODE[traj][(int) save_time][spec_i] = current_state[spec_i];
					}
					save_time += increment;
				}
			}

			// End of trajectory
			delete[] population_changes;
		}
	}
}
