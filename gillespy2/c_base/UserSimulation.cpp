#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include "model.h"
#include "ssa.h"
using namespace Gillespy;

bool seed_time = true;

//Default constants

const uint number_trajectories = 1;
const uint number_timesteps = 98;
const double end_time = 100000.0;
const double vol = 1.0;
int random_seed;
std :: string s_names[] = {"A", "B", "C", "X"};
uint populations[] = {300, 300, 300, 300};
std :: string r_names[] = {"j1", "j2"};
const double k1 = 1.0;
const double k2 = 1.0;

class PropensityFunction : public IPropensityFunction{
public:
  double evaluate(uint reaction_number, uint* state){
    switch(reaction_number){


        case 0:
            return k1*state[0]*state[3]/vol;

        

        case 1:
            return k2*state[1]*state[3]/vol;

        
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
model.reactions[0].species_change[0] = -1;
model.reactions[0].species_change[3] = 1.0;
model.reactions[1].species_change[1] = -1;
model.reactions[1].species_change[2] = 1;
model.reactions[1].species_change[3] = -1;
  //End reaction species changes
  model.update_affected_reactions();

 if(seed_time){
   random_seed = time(NULL);
 }
 
  IPropensityFunction *propFun = new PropensityFunction();
  Simulation simulation(&model, number_trajectories, number_timesteps, end_time, propFun, random_seed);
  ssa_direct(&simulation);
  std :: cout << simulation << std :: endl;
  delete propFun;
  return 0;
}
