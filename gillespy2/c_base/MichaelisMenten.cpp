#include <string>
#include <vector>
#include <iostream>
#include "model.h"
#include "ssa.h"
using namespace Gillespy;

const double rate1 = 0.0017;
const double rate2 = 0.5;
const double rate3 = 0.1;

const double volume = 1.0;


class PropensityFunction : public IPropensityFunction{
public:
  double evaluate(uint reaction_number, uint* state){
    switch(reaction_number){
    case 0:
      return rate1 * state[0] * state[1] / volume;
    case 1:
      return rate2 * state[2];
    case 2:
      return rate3 * state[2];
    default: //Error
      return -1;
    }
  }
};

int main(){
  std :: string species_names[] = {"A","B","C","D"};
  int initial_populations[] = {301, 120, 0, 0};
  std :: string reaction_names[] = {"r1","r2","r3"};
  std :: vector<std :: string> name_species, name_reactions;
  for(std :: string name : species_names){
    name_species.push_back(name);
  }
  for(std :: string name : reaction_names){
    name_reactions.push_back(name);
  }
  Model* model = build_model(name_species, name_reactions);
  for(int i = 0; i < model -> number_species; i++){
    model -> species[i].initial_population = initial_populations[i];
  }
  
  model -> reactions[0].species_change[0] = -1;
  model -> reactions[0].species_change[1] = -1;
  model -> reactions[0].species_change[2] = 1;

  model -> reactions[1].species_change[0] = 1;
  model -> reactions[1].species_change[1] = 1;
  model -> reactions[1].species_change[2] = -1;
  
  model -> reactions[2].species_change[1] = 1;
  model -> reactions[2].species_change[2] = -1;
  model -> reactions[2].species_change[3] = 1;

  for(uint i = 0; i < model -> number_reactions; i++){
    for(uint j = 0; j < model -> number_reactions; j++){
      model -> reactions[i].affected_reactions.push_back(j);
    }
  }
  
  IPropensityFunction *propFun = new PropensityFunction();
  Simulation* simulation = new Simulation(model, 7, 101, 100, propFun, 9001);
  ssa_direct(simulation);
  for(int i = 0; i < simulation -> number_timesteps; i++){
    std :: cout << simulation -> timeline[i] << ":\t";
    for(int j = 0; j < model -> number_species; j++){
      std :: cout << simulation -> trajectories[0][i][j] << ", ";
    }
    std :: cout << std :: endl;
  }
  delete simulation;
  delete propFun;
  free_model(model);
  return 0;
}
