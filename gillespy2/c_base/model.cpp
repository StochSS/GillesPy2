#include "model.h"

namespace Gillespy{
  
  Model* build_model(std :: vector<std :: string> name_species, std :: vector<std :: string> name_reactions){
    Model* model = new Model;
    model -> number_species = name_species.size();
    model -> number_reactions = name_reactions.size();
    model -> species = new Species[model -> number_species];
    for(uint i = 0; i < name_species.size(); i++){
      model -> species[i].id = i;
      model -> species[i].name = name_species[i];
    }
    model -> reactions = new Reaction[model -> number_reactions];
    for(uint i = 0; i < name_reactions.size(); i++){
      model -> reactions[i].name = name_reactions[i];
      model -> reactions[i].species_change = new int[model -> number_species];
      for(uint j = 0; j < model -> number_species; j++){
	model -> reactions[i].species_change[j] = 0;	
      }
      model -> reactions[i].affected_reactions = std :: vector<uint>();
      
    }
    return model;
  }

  void free_model(Model* model){
    if(model){
      if(model -> species){
	delete [] model -> species;
      }
      if(model -> reactions){
	for(uint i = 0; i < model -> number_reactions; i++){
	  if(model -> reactions[i].species_change){
	      delete [] model -> reactions[i].species_change;
	  }
	}
	delete [] model -> reactions;
      }
      delete model;
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
    free_model(this -> model);
    delete timeline;
    delete trajectories_1D;
    for(uint i = 0; i < number_trajectories; i++){
      delete trajectories[i];
    }
    delete trajectories;
  }
}
