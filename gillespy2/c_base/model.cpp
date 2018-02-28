#include "model.h"

namespace Gillespy{
  
  Model :: Model(std :: vector<std :: string> name_species, std :: vector<std :: string> name_reactions):
    number_species(name_species.size()),
    number_reactions(name_reactions.size())
  {
    species = std :: make_unique<Species[]>(number_species);
    for(uint i = 0; i < number_species; i++){
      species[i].id = i;
      species[i].name = name_species[i];
    }
    reactions = std :: make_unique<Reaction[]>(number_reactions);
    for(uint i = 0; i < number_reactions; i++){
      reactions[i].name = name_reactions[i];
      reactions[i].species_change = std :: make_unique<int[]>(number_species);
      for(uint j = 0; j < number_species; j++){
	reactions[i].species_change[j] = 0;	
      }
      reactions[i].affected_reactions = std :: vector<uint>();
    }
  }


  Simulation :: Simulation(Model* model, uint number_trajectories, uint number_timesteps, double end_time, IPropensityFunction* propensity_function, int random_seed) : model(model), end_time(end_time), random_seed(random_seed), number_timesteps(number_timesteps), number_trajectories(number_trajectories), propensity_function(propensity_function){
    timeline = new double[number_timesteps];
    double timestep_size = end_time/(number_timesteps-1);
    for(uint i = 0; i < number_timesteps; i++){
      timeline[i] = timestep_size * i;
    }
    uint trajectory_size = number_timesteps * (model -> number_species);
    trajectories_1D = new uint[number_trajectories * trajectory_size];
    trajectories = new uint**[number_trajectories];
    for(uint i = 0; i < number_trajectories; i++){
      trajectories[i] = new uint*[number_timesteps];
      for(uint j = 0; j < number_timesteps; j++){
	trajectories[i][j] = &(trajectories_1D[i * trajectory_size + j *  (model -> number_species)]);
      }
    }    
  }


  Simulation :: ~Simulation(){
    delete timeline;
    delete trajectories_1D;
    for(uint i = 0; i < number_trajectories; i++){
      delete trajectories[i];
    }
    delete trajectories;
  }
}
