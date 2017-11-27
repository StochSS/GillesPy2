#include "model.h"

namespace Gillespy{
  
  Model* build_model(int number_species, std :: vector<std :: string> name_species, int number_reactions, std :: vector<std :: string> name_reactions){
    Model* model = new Model;
    model -> species = new Species[number_species];
    for(int i = 0; i < number_species; i++){
      model -> species[i].id = i;
      model -> species[i].name = name_species[i];
    }
    model -> reactions = new Reaction[number_reactions];
    for(int i = 0; i < number_reactions; i++){
      model -> reactions[i].name = name_reactions[i];
      model -> reactions[i].species_change = new int[number_species];
      model -> reactions[i].affected_reactions = std :: vector<int>();
    }
    return model;
  }

  void free_model(Model* model){
    if(model){
      if(model -> species){
	delete model -> species;
      }
      if(model -> reactions){
	for(int i = 0; i < model -> number_reactions; i++){
	  if(model -> reactions[i].species_change){
	    delete model -> reactions[i].species_change;
	  }
	}
	delete model -> reactions;
      }
      delete model;
    }
  }

  Simulation :: Simulation(Model* model, int number_timesteps, int number_trajectories, double end_time, IPropensityFunction* propensity_function) : model(model), number_timesteps(number_timesteps), number_trajectories(number_trajectories), end_time(end_time), propensity_function(propensity_function){
    timeline = new double[number_timesteps];
    double timestep_size = end_time/(number_timesteps-1);
    for(int i = 0; i < number_timesteps; i++){
      timeline[i] = timestep_size * i;
    }
    int trajectory_size = number_timesteps * (model -> number_species);
    trajectories_1D = new int[number_trajectories * trajectory_size];
    trajectories = new int**[number_trajectories];
    for(int i = 0; i < number_trajectories; i++){
      trajectories[i] = new int*[number_timesteps];
      for(int j = 0; j < number_timesteps; j++){
	trajectories[i][j] = &(trajectories_1D[i * trajectory_size + j *  (model -> number_species)]);
      }
    }    
  }


  Simulation :: ~Simulation(){
    free_model(this -> model);
    delete timeline;
    delete trajectories_1D;
    for(int i = 0; i < number_trajectories; i++){
      delete trajectories[i];
    }
    delete trajectories;
  }
}
