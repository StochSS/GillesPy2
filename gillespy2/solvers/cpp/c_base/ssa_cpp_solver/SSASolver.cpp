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

#include <cmath>
#include <random>
#include "SSASolver.h"

namespace Gillespy
{
    volatile bool interrupted = false;

    GPY_INTERRUPT_HANDLER(signal_handler, {
        interrupted = true;
    })

    void ssa_direct(Simulation<unsigned int> *simulation)
    {
        GPY_INTERRUPT_INSTALL_HANDLER(signal_handler);

        if (simulation)
        {
            std::mt19937_64 rng(simulation->random_seed);

            //Number of bytes for copying states
            unsigned int state_size = sizeof(int) * ((simulation->model)->number_species);

            //Calculated propensity values for current state
            double *propensity_values = new double[(simulation->model)->number_reactions];

            //Simulate for each trajectory
            for (unsigned int trajectory_number = 0; trajectory_number < simulation->number_trajectories; trajectory_number++)
            {
                if (interrupted)
                {
                    break;
                }

                //Set up current state from initial state
                unsigned int entry_count = 0;
                simulation->current_time = 0;
                simulation->reset_output_buffer(trajectory_number);
                simulation->output_buffer_range(std::cout);

                //calculate initial propensities
                for (unsigned int reaction_number = 0; reaction_number < ((simulation->model)->number_reactions); reaction_number++)
                {
                    propensity_values[reaction_number] = Reaction::propensity(reaction_number, simulation->current_state);
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
                        //Quit simulating this trajectory
                        break;
                    }

                    //Reaction will fire, determine which one
                    double cumulative_sum = rng() * propensity_sum / rng.max();
                    simulation->current_time += -log(rng() * 1.0 / rng.max()) / propensity_sum;

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
                                simulation->current_state[species_number] += reaction.species_change[species_number];
                            }

                            //Recalculate needed propensities
                            for (unsigned int &affected_reaction : reaction.affected_reactions)
                            {
                                propensity_values[affected_reaction] = Reaction::propensity(affected_reaction, simulation->current_state);
                            }

                            break;
                        }
                    }

                    // Output timestep results
                    while (entry_count < simulation->number_timesteps && simulation->timeline[entry_count] <= simulation->current_time)
                    {
                        if (interrupted)
                        {
                            break;
                        }

                        simulation->output_buffer_range(std::cout, entry_count++);
                    }
                }

                
                // Copy final state for rest of entries
                if (entry_count < simulation->number_timesteps)
                {
                    simulation->current_time = simulation->timeline[simulation->number_timesteps - 1];
                    simulation->output_buffer_range(std::cout, simulation->number_timesteps - 1);
                }
            }

            delete propensity_values;
        }
    }
}
