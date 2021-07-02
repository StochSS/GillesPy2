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

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include <time.h>
#include <math.h>

#include "TauLeapingSolver.h"
#include "template.h"
#include "arg_parser.h"

using namespace Gillespy;

bool seed_time = true;

int random_seed = 0;
unsigned int number_trajectories = 0;
unsigned int number_timesteps = 0;

double end_time = 0;
double tau_tol = 0.03;

class PropensityFunction : public IPropensityFunction
{
public:
	double TauEvaluate(unsigned int reaction_number, const std::vector<int> &S)
	{
		return map_propensity(reaction_number, S);
	}

	double evaluate(unsigned int reaction_number, unsigned int *state)
	{
		return 1.0;
	}

	double ODEEvaluate(int reaction_number, const std::vector<double> &S)
	{
		return 1.0;
	}
};

int main(int argc, char* argv[]){
	ArgParser parser = ArgParser(argc, argv);

	random_seed = parser.seed;

	if (random_seed != -1)
	{
		seed_time = false;
	}

	end_time = parser.end;
	number_trajectories = parser.trajectories;
	number_timesteps = parser.timesteps;
	tau_tol = parser.tau_tol;

	Model<unsigned int> model(species_names, species_populations, reaction_names);
	add_reactions(model);

	if(seed_time) {
		random_seed = time(NULL);
	}

	IPropensityFunction *propFun = new PropensityFunction();
	Simulation<unsigned int> simulation;

	simulation.model = &model;
	simulation.end_time = end_time;
	simulation.random_seed = random_seed;
	simulation.number_timesteps = number_timesteps;
	simulation.number_trajectories = number_trajectories;
	simulation.propensity_function = propFun;

	init_simulation(&model, simulation);
	tau_leaper(&simulation, tau_tol);

	simulation.output_results_buffer(std :: cout);

	delete propFun;
	return 0;
}
