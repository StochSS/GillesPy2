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
#include <iostream>
#include <sstream>
#include <time.h>
#include <math.h>
#include "TauHybridSolver.h"
#include "template.h"
#include "hybrid_template.h"
#include "arg_parser.h"
using namespace Gillespy;

//Default values, replaced with command line args
unsigned int number_trajectories = 0;
unsigned int number_timesteps = 0;
int random_seed = 0;
double end_time = 100.0;
bool seed_time = true;
double increment = 0;
double tau_tol = 0.03;

int main(int argc, char* argv[])
{
    ArgParser parser(argc, argv);
    random_seed = parser.seed;
    seed_time = (random_seed == -1);

    end_time = parser.end;
    number_trajectories = parser.trajectories;
    number_timesteps = parser.timesteps;
    tau_tol = parser.tau_tol;

    Reaction::load_parameters();
    Model<double> model(species_names, species_populations, reaction_names);
    add_reactions(model);

    if(seed_time){
        random_seed = time(nullptr) % GPY_PID_GET();
    }
    //Simulation INIT
    TauHybrid::HybridSimulation simulation(model);
    simulation.model = &model;
    simulation.end_time = end_time;
    simulation.random_seed = random_seed;
    simulation.number_timesteps = number_timesteps;
    simulation.number_trajectories = number_trajectories;
    simulation.output_interval = parser.output_interval;
    init_simulation(&model, simulation);
    Gillespy::TauHybrid::map_species_modes(simulation.species_state);
    Gillespy::TauHybrid::map_rate_rules(simulation.species_state);

    std::vector<Gillespy::TauHybrid::Event> events;
    Gillespy::TauHybrid::Event::use_events(events);
    Logger logger;
    if (parser.verbose)
        logger.set_log_level(LogLevel::INFO);

    // Configure solver options from command-line arguments (or defaults)
    SolverConfiguration config = {
        parser.rtol,
        parser.atol,
        parser.max_step,

    };

    TauHybrid::TauHybridCSolver(&simulation, events, logger, tau_tol, config, parser.use_root_finding);
    simulation.output_buffer_final(std::cout);
    return simulation.get_status();
}
