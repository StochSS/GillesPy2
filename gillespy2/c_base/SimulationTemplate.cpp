#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include "model.h"
#include "ssa.h"
using namespace Gillespy;

bool seed_time = true;

//Default constants
__DEFINE_CONSTANTS__

class PropensityFunction : public IPropensityFunction{
public:
  double evaluate(uint reaction_number, uint* state){
    switch(reaction_number){
__DEFINE_PROPENSITY__

    default: //Error
      return -1;
    }
  }
};

int main(){
  std :: vector<std :: string> species_names(s_names, s_names + sizeof(s_names)/sizeof(std :: string));
  std :: vector<uint> species_populations(populations, populations + sizeof(populations)/sizeof(populations[0]));
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
  Simulation simulation(&model, number_trajectories, number_timesteps, end_time, propFun, random_seed);
  ssa_direct(&simulation);
  //std :: cout << simulation << std :: endl;
  simulation.output_results_buffer(std :: cout);
  delete propFun;
  return 0;
}
