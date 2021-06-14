#include <cmath>
#include <random>
#include <csignal>

#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include "SSASolver.h"

namespace Gillespy {
	volatile bool interrupted = false;

	#ifdef _WIN32
	BOOL WINAPI eventHandler(DWORD CtrlType) {
		interrupted = true;
		return TRUE;
	}
	#endif

	void signalHandler(int signum) {
		interrupted = true;
	}

	void ssa_direct(Simulation<unsigned int> *simulation) {
		#ifdef _WIN32
		SetConsoleCtrlHandler(eventHandler, TRUE);
		#else
		signal(SIGINT, signalHandler);
		#endif

		Model *model = simulation->model;

		if (!simulation) {
			return;
		}

		std::mt19937_64 rng(simulation->random_seed);

		unsigned int state_size = sizeof(int) * model->number_species;
		unsigned int *current_state = new unsigned int[model->number_species];

		double *propensity_values = new double[model->number_reactions];

		// Copy initial populations into the trajectories array.
		for (unsigned int species_index = 0; species_index < model->number_reactions; species_index++) {
			simulation->trajectories[0][0][species_index] = model->species[species_index].initial_population;
		}

		// Begin simulating each trajectory.
		for (unsigned int trajectory_index = 0; trajectory_index < simulation->number_trajectories; trajectory_index++) {
			if (interrupted) {
				break;
			}

			// Grab a references to the current trajectory.
			unsigned int **trajectory = simulation->trajectories[trajectory_index];


			// If needed, copy the initial state.
			if (trajectory_index > 0) {
				memcpy(trajectory[0], simulation->trajectories[0][0], state_size);
			}

			// Setup the current state.
			memcpy(current_state, trajectory[0], state_size);
			simulation->current_time = 0;

			unsigned int entry_count = 1;

			// Calculate initial propensities.
			for (unsigned int reaction_index = 0; reaction_index < model->number_reactions; reaction_index++) {
				propensity_values[reaction_index] = simulation->propensity_function->evaluate(reaction_index, current_state);
			}

			double propensity_sum;

			while (simulation->current_time < simulation->end_time) {
				if (interrupted) {
					break;
				}

				propensity_sum = 0;
				for (unsigned int reaction_index = 0; reaction_index < model->number_reactions; reaction_index++) {
					propensity_sum += propensity_values[reaction_index];
				}

				// If there are no more reactions, copy the last state into the trajectories.
				if (propensity_sum <= 0) {
					for (unsigned int i = entry_count; i < simulation->number_timesteps; i++) {
						memcpy(trajectory[i], current_state, state_size);
					}

					break;
				}

				// Since a reaction must fire, determine which.
				double cumulative_sum = rng() * propensity_sum / rng.max();
				simulation->current_time += -log(rng() * 1.0 / rng.max()) / propensity_sum;

				// Copy the current state into the trajectory.
				while (entry_count < simulation->number_timesteps && simulation->timeline[entry_count] <= simulation->current_time) {
					if (interrupted) {
						break;
					}

					memcpy(trajectory[entry_count], current_state, state_size);
					entry_count++;
				}

				for (unsigned int potential_reaction = 0; potential_reaction < model->number_reactions; potential_reaction++) {
					cumulative_sum -= propensity_values[potential_reaction];

					// If cumulative_sum <= 0, it means that the reaction at potential_reaction has been fired.
					if (cumulative_sum <= 0 && propensity_values[potential_reaction] > 0) {
						// Update the state of the current reaction.
						Reaction &reaction = model->reactions[potential_reaction];

						for (unsigned int species_index = 0; species_index < model->number_species; species_index++) {
							current_state[species_index] += reaction.species_change[species_index];
						}

						// Recalculate any subsequent dependent propensities.
						for (unsigned int &affected_reaction : reaction.affected_reactions) {
							propensity_values[affected_reaction] = simulation->propensity_function->evaluate(affected_reaction, current_state);
						}

						break;
					}
				}
			}
		}

		// Cleanup and have a nice day!
		delete[] propensity_values;
		delete[] current_state;
	}
}
