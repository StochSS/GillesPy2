#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include <time.h>
#include <math.h>

#include "TauLeapingSolver.h"
#include "template.h"

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

int main(int argc, char *argv[])
{
	//Parse command line arguments
	// TODO: NEEDS REPLACEMENT
	std::string arg;
	for (int i = 1; i < argc - 1; i++)
	{
		arg = argv[i];
		if (argc > i + 1 && arg.size() > 1 && arg[0] == '-')
		{
			std::stringstream arg_stream(argv[i + 1]);
			switch (arg[1])
			{
			case 's':
				arg_stream >> random_seed;
				seed_time = false;
				break;
			case 'e':
				arg_stream >> end_time;
				break;
			case 'i':
				map_variable_populations(arg_stream);
				break;
			case 'p':
				map_variable_parameters(arg_stream);
				break;
			case 't':
				if (arg[2] == 'r')
				{
					arg_stream >> number_trajectories;
				}
				else if (arg[2] == 'i')
				{
					arg_stream >> number_timesteps;
				}

				else if (arg[2] == 'a')
				{ // '-tau_tol'
					arg_stream >> tau_tol;
				}
				break;
			}
		}
	}

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

	tau_leaper(&simulation, tau_tol);
	simulation.output_results_buffer(std::cout);

	delete propensity_function;
	return 0;
}
