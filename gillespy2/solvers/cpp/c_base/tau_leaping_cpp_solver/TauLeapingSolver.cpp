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

#include <map>
#include <set>
#include <cmath>
#include <vector>
#include <cstring>
#include <string>
#include <random>
#include <csignal>
#include <iostream>
#include <functional>
#include <algorithm>

#include "TauLeapingSolver.h"
#include "tau.h"

namespace Gillespy
{
    static volatile bool interrupted = false;
    std::mt19937_64 generator;

    GPY_INTERRUPT_HANDLER(signal_handler, {
        interrupted = true;
    })

    std::pair<std::map<std::string, int>, double> get_reactions(
        const Gillespy::Model<unsigned int> *model,
        const std::vector<double> &propensity_values,
        double tau_step,
        double current_time,
        double save_time)
    {

        /*
        * Helper Function to get reactions fired from t to t+tau. Effects two values:
        *rxn_count - dict with key=Reaction channel value=number of times fired
        *curr_time - float representing current time
        */

        if (current_time + tau_step > save_time)
        {
            tau_step = save_time - current_time;
        }

        std::map<std::string, int> rxn_count; // map of how many times reaction is fired
        std::pair<std::map<std::string, int>, double> values; // value pair to be returned, map of times {map of times reaction fired, current time}

        for (int i = 0; i < model->number_reactions; i++)
        {
            std::poisson_distribution<int> poisson(propensity_values[i] * tau_step);
            rxn_count[model->reactions[i].name] = poisson(generator);
        }

        current_time = current_time + tau_step;

        values.first = rxn_count;
        values.second = current_time;
        return values;
    }

    void tau_leaper(Gillespy::Simulation<unsigned int> *simulation, const double tau_tol)
    {
        GPY_INTERRUPT_INSTALL_HANDLER(signal_handler);

        if (!simulation)
        {
            return;
        }

        //Initialize your tau args
        TauArgs<unsigned int> tau_args = initialize(*(simulation->model), tau_tol);

        double increment = simulation->timeline[1] - simulation->timeline[0];

        // Instantiate the RNG.
        generator = std::mt19937_64(simulation->random_seed);

        //Initialize current_state variables, propensity_values
        std::vector<int> current_state(simulation->model->number_species);
        std::vector<double> propensity_values(simulation->model->number_reactions);

        //Simulate for each trajectory
        for (unsigned int trajectory_number = 0; trajectory_number < simulation->number_trajectories; trajectory_number++)
        {
            if (interrupted)
            {
                break;
            }

            simulation->reset_output_buffer(trajectory_number);
            //copy initial state for each trajectory
            for (unsigned int species_number = 0; species_number < (simulation->model->number_species); species_number++)
            {
                current_state[species_number] = simulation->model->species[species_number].initial_population;
            }
            std::copy(current_state.begin(), current_state.end(), simulation->current_state);

            //Initialize simulation variables
            simulation->current_time = 0;
            unsigned int entry_count = 0;

            //Propensity sum initialization, to be added to later.
            double propensity_sum;

            //Start save time
            double save_time = 0;

            //Variable to keep track of rejected steps, debug
            int steps_rejected = 0;

            //Initialize tau_step, will be assigned using tau::select()
            double tau_step;
            std::vector<int> prev_curr_state(current_state.begin(), current_state.end());

            // Each save step
            while (entry_count < simulation->number_timesteps)
            { // while less than end_time? Could be incorrect
                if (interrupted)
                {
                    break;
                }

                while (simulation->current_time < save_time)
                {
                    if (interrupted)
                    {
                        break;
                    }

                    //calculate propensities for each step
                    for (unsigned int reaction_number = 0; reaction_number < simulation->model->number_reactions; reaction_number++)
                    {
                        propensity_values[reaction_number] = Reaction::propensity(reaction_number, current_state.data());
                    }

                    tau_step = select(*(simulation->model), tau_args, tau_tol, simulation->current_time, save_time, propensity_values, current_state);
                    prev_curr_state = current_state;
                    double prev_curr_time = simulation->current_time;
                    int loop_cnt = 0;

                    while (true)
                    {
                        loop_cnt += 1;

                        if (loop_cnt > 100)
                        {
                            throw std::runtime_error("Loop count exceeded 100, error");
                        }

                        std::map<std::string, int> rxn_count;
                        std::pair<std::map<std::string, int>, double> values;

                        values = get_reactions(simulation->model, propensity_values, tau_step, simulation->current_time, save_time);

                        rxn_count = values.first;
                        simulation->current_time = values.second;

                        std::map<int, bool> species_modified;

                        for (int i = 0; i < simulation->model->number_reactions; i++)
                        {
                            if (rxn_count[simulation->model->reactions[i].name] > 0)
                            {
                                for (auto const &spec : tau_args.reactions_reactants[i])
                                {
                                    species_modified[spec] = true;
                                    //+= for both reactants and products because current_state is represented with negative number changes for reactants, and positive for products.
                                    current_state[spec] += simulation->model->reactions[i].species_change[spec] * rxn_count[simulation->model->reactions[i].name];
                                }
                            }

                            for (auto const &spec : tau_args.products[i])
                            {
                                species_modified[spec] = true;
                                current_state[spec] += simulation->model->reactions[i].species_change[spec] * rxn_count[simulation->model->reactions[i].name];
                            }
                        }


                        bool neg_state = false;
                        for (auto const &x : species_modified)
                        {
                            if (current_state[x.first] < 0)
                            {
                                neg_state = true;
                            }
                        }

                        if (neg_state)
                        {
                            current_state = prev_curr_state;
                            simulation->current_time = prev_curr_time;
                            tau_step /= 2;
                        }

                        else
                        {
                            break; // out of while true
                        }
                    }
                }
                // Copy internal vector into simulation state array
                std::copy(current_state.begin(), current_state.end(), simulation->current_state);
                simulation->output_buffer_range(std::cout);

                save_time += increment;
                entry_count += 1;
            }
        }
    }
}
