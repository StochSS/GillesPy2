#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <time.h>
#include <math.h>
#include "model.h"
#include "ssa.h"
using namespace Gillespy;

//Default values, replaced with command line args
unsigned int number_trajectories = 0;
unsigned int number_timesteps = 0;
int random_seed = 0;
double end_time = 0;
bool seed_time = true;

//Default constants/variables
__DEFINE_VARIABLES__

class PropensityFunction : public IPropensityFunction{
public:
  double evaluate(unsigned int reaction_number, unsigned int* S){
    switch(reaction_number){
__DEFINE_PROPENSITY__

    default: //Error
      return -1;
    }
  }
  double TauEvaluate(unsigned int reaction_number,  const std::vector<int> &S){return 1.0;}
  double ODEEvaluate(int reaction_number, const std::vector <double> &S){return 1.0;}
};

int main(int argc, char* argv[]){
 //Parse command line arguments
 std :: string arg;
 for(int i = 1; i < argc - 1; i++){
   arg = argv[i];
   if(argc > i+1 && arg.size() > 1 && arg[0] == '-'){
     std :: stringstream arg_stream(argv[i+1]);
     switch(arg[1]){
     case 'i':
       for(int j = 0; j < int(sizeof(populations)); j++){
       arg_stream >> populations[j];
       }
       break;
     case 'p':
__DEFINE_PARAMETER_UPDATES__
       break;
     case 's':
       arg_stream >> random_seed;
       seed_time = false;
       break;
     case 'e':
       arg_stream >> end_time;
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

  std :: vector<std :: string> species_names(s_names, s_names + sizeof(s_names)/sizeof(std :: string));
  std :: vector<unsigned int> species_populations(populations, populations + sizeof(populations)/sizeof(populations[0]));
  std :: vector<std :: string> reaction_names(r_names, r_names + sizeof(r_names)/sizeof(std :: string));
  
  Model model(species_names, species_populations, reaction_names);

  //Begin reaction species changes
__DEFINE_REACTIONS_
  //End reaction species changes
  model.update_affected_reactions();
 
 if(seed_time){
   random_seed = time(NULL);
 }
  IPropensityFunction *propFun = new PropensityFunction();
  // Simulation INIT
  Simulation simulation;
  Model* modelptr;
  modelptr = &model;
  simulation.model = modelptr;
  simulation.end_time = end_time;
  simulation.random_seed = random_seed;
  simulation.number_timesteps = number_timesteps;
  simulation.number_trajectories = number_trajectories;
  simulation.propensity_function = propFun;
  simulationSSAINIT(&model, simulation);
  // Perform SSA  //
  ssa_direct(&simulation);

  simulation.output_results_buffer(std :: cout);
  delete propFun;
  return 0;
}
