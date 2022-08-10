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

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include <time.h>
#include <math.h>

#include "SSASolver.h"
#include "template.h"
#include "arg_parser.h"

using namespace Gillespy;

bool seed_time = true;

int random_seed = 0;
unsigned int number_trajectories = 0;
unsigned int number_timesteps = 0;

double end_time = 0;

int main(int argc, char *argv[])
{
    //Parse command line arguments
    ArgParser parser = ArgParser(argc, argv);
    
    random_seed = parser.seed;
    if (random_seed != -1)
    {
        seed_time = false;
    }
    
    end_time = parser.end;
    number_trajectories = parser.trajectories;
    number_timesteps = parser.timesteps;

    Reaction::load_parameters();
    Model<unsigned int> model(species_names, species_populations, reaction_names);
    add_reactions(model);
    
    if (seed_time)    
    {
        random_seed = time(nullptr) % GPY_PID_GET();
    }

    Simulation<unsigned int> simulation;
    simulation.output_interval = parser.output_interval;
    simulation.model = &model;
    simulation.end_time = end_time;
    simulation.random_seed = random_seed;
    simulation.number_timesteps = number_timesteps;
    simulation.number_trajectories = number_trajectories;
    
    init_simulation(&model, simulation);
    ssa_direct(&simulation);
    simulation.output_buffer_final(std::cout);

    return simulation.get_status();
}
