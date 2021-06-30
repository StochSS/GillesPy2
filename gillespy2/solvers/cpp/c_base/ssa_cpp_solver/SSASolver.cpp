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

#include <cmath>
#include <random>
#include <csignal>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include "SSASolver.h"

namespace Gillespy
{
	volatile bool interrupted = false;

#ifdef _WIN32
	BOOL WINAPI eventHandler(DWORD CtrlType)
	{
		interrupted = true;
		return TRUE;
	}
#endif

	void signalHandler(int signum)
	{
		interrupted = true;
	}

	void ssa_direct(Simulation<unsigned int> *simulation)
	{
#ifdef _WIN32
		SetConsoleCtrlHandler(eventHandler, TRUE);
#else
		signal(SIGINT, signalHandler);
#endif

		if (simulation)
		{
			std::mt19937_64 rng(simulation->random_seed);

			//Number of bytes for copying states
			unsigned int state_size = sizeof(int) * ((simulation->model)->number_species);

			//Current state
			unsigned int *current_state = new unsigned int[(simulation->model)->number_species];

			//Calculated propensity values for current state
			double *propensity_values = new double[(simulation->model)->number_reactions];

			//copy initial state for each trajectory
			for (unsigned int species_number = 0; species_number < ((simulation->model)->number_species); species_number++)
			{
				simulation->trajectories[0][0][species_number] = (simulation->model)->species[species_number].initial_population;
			}

			//Simulate for each trajectory
			for (unsigned int trajectory_number = 0; trajectory_number < simulation->number_trajectories; trajectory_number++)
			{
				if (interrupted)
				{
					break;
				}

				//Get simpler reference to memory space for this trajectory
				unsigned int **trajectory = simulation->trajectories[trajectory_number];

				//Copy initial state as needed
				if (trajectory_number > 0)
				{
					memcpy(trajectory[0], simulation->trajectories[0][0], state_size);
				}

				//Set up current state from initial state
				memcpy(current_state, trajectory[0], state_size);
				simulation->current_time = 0;
				unsigned int entry_count = 1;

				//calculate initial propensities
				for (unsigned int reaction_number = 0; reaction_number < ((simulation->model)->number_reactions); reaction_number++)
				{
					propensity_values[reaction_number] = (simulation->propensity_function)->evaluate(reaction_number, current_state);
				}

				double propensity_sum;
				while (simulation->current_time < (simulation->end_time))
				{
					if (interrupted)
					{
						break;
					}

					//Sum propensities
					propensity_sum = 0;
					for (unsigned int reaction_number = 0; reaction_number < ((simulation->model)->number_reactions); reaction_number++)
					{
						propensity_sum += propensity_values[reaction_number];
					}

					//No more reactions
					if (propensity_sum <= 0)
					{
						//Copy all of last changed state for rest of entries
						for (unsigned int i = entry_count; i < simulation->number_timesteps; i++)
						{
							memcpy(trajectory[i], current_state, state_size);
						}

						//Quit simulating this trajectory
						break;
					}

					//Reaction will fire, determine which one
					double cumulative_sum = rng() * propensity_sum / rng.max();
					simulation->current_time += -log(rng() * 1.0 / rng.max()) / propensity_sum;

					//Copy current state to passed timesteps
					while (entry_count < simulation->number_timesteps && (simulation->timeline[entry_count]) <= simulation->current_time)
					{
						if (interrupted)
						{
							break;
						}

						memcpy(trajectory[entry_count], current_state, state_size);
						entry_count++;
					}

					for (unsigned int potential_reaction = 0; potential_reaction < ((simulation->model)->number_reactions); potential_reaction++)
					{
						cumulative_sum -= propensity_values[potential_reaction];
						//This reaction fired

						if (cumulative_sum <= 0 && propensity_values[potential_reaction] > 0)
						{
							//Update current state
							Reaction &reaction = ((simulation->model)->reactions[potential_reaction]);
							for (unsigned int species_number = 0; species_number < ((simulation->model)->number_species); species_number++)
							{
								current_state[species_number] += reaction.species_change[species_number];
							}

							//Recalculate needed propensities
							for (unsigned int &affected_reaction : reaction.affected_reactions)
							{
								propensity_values[affected_reaction] = (simulation->propensity_function)->evaluate(affected_reaction, current_state);
							}

							break;
						}
					}
				}
			}

			delete propensity_values;
			delete current_state;
		}
	}
}
