#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include <time.h>
#include <math.h>

#include "ODESolver.h"
#include "template.h"

using namespace Gillespy;

bool seed_time = true;

int random_seed = 0;
unsigned int number_trajectories = 0;
unsigned int number_timesteps = 0;

double end_time = 100.0;
double increment = 0;

class PropensityFunction : public IPropensityFunction {
public:
	double evaluate(unsigned int reaction_number, unsigned int *S) {
		return 1.0;
	}

	double TauEvaluate(unsigned int reaction_number, const std::vector<int> &S) {
		return 1.0;
	}

	double ODEEvaluate(int reaction_number, const std::vector<double> &S) {
		return map_ode_propensity(reaction_number, S);
	}
};

int main(int argc, char *argv[]) {
//Parse command line arguments
// TODO: NEEDS REPLACEMENT.
 std :: string arg;
 for(int i = 1; i < argc - 1; i++){
   arg = argv[i];
    if(argc > i+1 && arg.size() > 1 && arg[0] == '-'){
        std :: stringstream arg_stream(argv[i+1]);
        switch(arg[1]){
        case 's':
            arg_stream >> random_seed;
            seed_time = false;
            break;
        case 'e':
            arg_stream >> end_time;
            break;
        case 'i':
            if (arg[3] == 'c') {
                arg_stream >> increment;
            }
            else if (arg[3] == 'i') {
                map_variable_populations(arg_stream);
            }
            break;
        case 'p':
            map_variable_parameters(arg_stream);
            break;
        case 't':
            if(arg[2] == 'r'){
                arg_stream >> number_trajectories;
            }else if(arg[2] == 'i'){
                arg_stream >> number_timesteps;
            }
            break;
        }
   }
 }

	Model model(species_names, species_populations, reaction_names);
	add_reactions(model);

	if (seed_time) {
		random_seed = time(NULL);
	}

	IPropensityFunction *propensity_function = new PropensityFunction();
	Simulation<double> simulation;

	simulation.model = &model;
	simulation.end_time = end_time;
	simulation.random_seed = random_seed;
	simulation.number_timesteps = number_timesteps;
	simulation.number_trajectories = number_trajectories;
	simulation.propensity_function = propensity_function;

	init_simulation(&model, simulation);

	ODESolver(&simulation, increment);
	simulation.output_results_buffer(std::cout);

	delete propensity_function;

	return 0;
 }