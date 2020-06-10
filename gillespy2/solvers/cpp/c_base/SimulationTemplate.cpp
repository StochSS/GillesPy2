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

//Default constants
__DEFINE_CONSTANTS__

class PropensityFunction : public IPropensityFunction{
public:
  double evaluate(unsigned int reaction_number, unsigned int* S){
    switch(reaction_number){
__DEFINE_PROPENSITY__

    default: //Error
      return -1;
    }
  }
};

int main(int argc, char* argv[]){
  std :: vector<std :: string> species_names(s_names, s_names + sizeof(s_names)/sizeof(std :: string));
  std :: vector<unsigned int> species_populations(populations, populations + sizeof(populations)/sizeof(populations[0]));
  std :: vector<std :: string> reaction_names(r_names, r_names + sizeof(r_names)/sizeof(std :: string));
  
  Model model(species_names, species_populations, reaction_names);

  //Begin reaction species changes
__DEFINE_REACTIONS_
  //End reaction species changes
  model.update_affected_reactions();
 
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

 if(seed_time){
   random_seed = time(NULL);
 }
  IPropensityFunction *propFun = new PropensityFunction();
  Simulation simulation(&model, number_trajectories, number_timesteps, end_time, propFun, random_seed);
  ssa_direct(&simulation);
  //std :: cout << simulation << std :: endl;
  simulation.output_results_buffer(std :: cout);
  delete propFun;
  return 0;
}
