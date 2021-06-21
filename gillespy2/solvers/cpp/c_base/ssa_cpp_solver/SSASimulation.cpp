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

class PropensityFunction : public IPropensityFunction
{
public:
	double evaluate(unsigned int reaction_number, unsigned int *S)
	{
		return map_propensity(reaction_number, S);
	}

	double TauEvaluate(unsigned int reaction_number, const std::vector<int> &S)
	{
		return 1.0;
	}

	double ODEEvaluate(int reaction_number, const std::vector <double> &S)
	{
		return 1.0;
	}
};

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
	
	Model model(species_names, species_populations, reaction_names);
	add_reactions(model);
	
	if (seed_time)    
	{
		random_seed = time(NULL);
	}
	
	IPropensityFunction *propensity_function = new PropensityFunction();
	Simulation<unsigned int> simulation;
	
	simulation.model = &model;
	simulation.end_time = end_time;
	simulation.random_seed = random_seed;
	simulation.number_timesteps = number_timesteps;
	simulation.number_trajectories = number_trajectories;
	simulation.propensity_function = propensity_function;
	
	init_simulation(&model, simulation);
	
	ssa_direct(&simulation);
	simulation.output_results_buffer(std::cout);
	
	delete propensity_function;
	return 0;
}
