#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <time.h>
#include <math.h>
#include "ODESolver.h"
#include "template.h" // User-defined simulation parameters

using namespace Gillespy;

//Default values, replaced with command line args
unsigned int number_trajectories = 0;
unsigned int number_timesteps = 0;
int random_seed = 0;
double end_time = 100.0;
bool seed_time = true;
double increment = 0;

class PropensityFunction : public IPropensityFunction{
public:
  double ODEEvaluate(int reaction_number, const std::vector <double> &S){
    return map_ode_propensity(reaction_number, S);
  }
  double TauEvaluate(unsigned int reaction_number, const std::vector<int> &S){return 1.0;}
  double evaluate(unsigned int reaction_number, unsigned int* S){return 1.0;}
};

int main(int argc, char* argv[]){
  //Parse command line arguments
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
  // Reactions are defined via macros in the user-defined template.
  // Effectively translates from a 2d array to the model's internal data structure.
  add_reactions(model);

 if(seed_time){
   random_seed = time(NULL);
 }
  IPropensityFunction *propFun = new PropensityFunction();
   //Simulation INIT
  Simulation simulation;
  Model* modelptr;
  modelptr = &model;
  simulation.ISODE=1;
  simulation.model = modelptr;
  simulation.end_time = end_time;
  simulation.random_seed = random_seed;
  simulation.number_timesteps = number_timesteps;
  simulation.number_trajectories = number_trajectories;
  simulation.propensity_function = propFun;
  simulationODEINIT(&model, simulation);
  // Perform ODE  //
  ODESolver(&simulation,increment);
  simulation.output_results_buffer(std :: cout);
  delete propFun;
  return 0;
}