#include "model.h"

namespace Gillespy{
  
  Model :: Model(std :: vector<std :: string> species_names, std :: vector<uint> species_populations, std :: vector<std :: string> reaction_names):
    number_species(species_names.size()),
    number_reactions(reaction_names.size())
  {
    species = std :: make_unique<Species[]>(number_species);
    for(uint i = 0; i < number_species; i++){
      species[i].id = i;
      species[i].initial_population = species_populations[i];
      species[i].name = species_names[i];
    }
    reactions = std :: make_unique<Reaction[]>(number_reactions);
    for(uint i = 0; i < number_reactions; i++){
      reactions[i].name = reaction_names[i];
      reactions[i].species_change = std :: make_unique<int[]>(number_species);
      for(uint j = 0; j < number_species; j++){
	reactions[i].species_change[j] = 0;	
      }
      reactions[i].affected_reactions = std :: vector<uint>();
    }
  }

  void Model :: update_affected_reactions(){
    //Clear affected_reactions for each reaction
    for(uint i = 0; i < number_reactions; i++){
      reactions[i].affected_reactions.clear();
    }   
    //Check all reactions for common species changes -> affected reactions
    for(uint r1 = 0; r1 < number_reactions; r1++){
      reactions[r1].affected_reactions.push_back(r1);
      for(uint r2 = r1 + 1; r2 < number_reactions; r2++){
	for(uint s = 0; s < number_species; s++){
	  if(reactions[r1].species_change[s] != 0 and reactions[r2].species_change[s] != 0){
	    reactions[r1].affected_reactions.push_back(r2);
	    reactions[r2].affected_reactions.push_back(r1);
	  }
	}
      }
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

  
  std :: ostream& operator<<(std :: ostream& os, const Simulation& simulation){
    for(uint i = 0; i < simulation.number_timesteps; i++){
      os << simulation.timeline[i] << " ";
      for(uint trajectory = 0; trajectory < simulation.number_trajectories; trajectory++){
	for(uint j = 0; j < simulation.model -> number_species; j++){
	  os << simulation.trajectories[trajectory][i][j] <<  " ";
	}
      }
      os << "\n";
    }
    return os;
  }
}
