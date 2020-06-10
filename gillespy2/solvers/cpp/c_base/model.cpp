#include "model.h"

namespace Gillespy{
  
  Model :: Model(std :: vector<std :: string> species_names, std :: vector<unsigned int> species_populations, std :: vector<std :: string> reaction_names):
    number_species(species_names.size()),
    number_reactions(reaction_names.size())
  {
    species = std :: make_unique<Species[]>(number_species);
    for(unsigned int i = 0; i < number_species; i++){
      species[i].id = i;
      species[i].initial_population = species_populations[i];
      species[i].name = species_names[i];
    }
    reactions = std :: make_unique<Reaction[]>(number_reactions);
    for(unsigned int i = 0; i < number_reactions; i++){
      reactions[i].name = reaction_names[i];
      reactions[i].species_change = std :: make_unique<int[]>(number_species);
      for(unsigned int j = 0; j < number_species; j++){
	reactions[i].species_change[j] = 0;	
      }
      reactions[i].affected_reactions = std :: vector<unsigned int>();
    }
  }

  void Model :: update_affected_reactions(){
    //Clear affected_reactions for each reaction
    for(unsigned int i = 0; i < number_reactions; i++){
      reactions[i].affected_reactions.clear();
    }   
    //Check all reactions for common species changes -> affected reactions
    for(unsigned int r1 = 0; r1 < number_reactions; r1++){
      for(unsigned int r2 = 0; r2 < number_reactions; r2++){
	    for(unsigned int s = 0; s < number_species; s++){
	        if(reactions[r2].species_change[s] != 0){
	            reactions[r1].affected_reactions.push_back(r2);
	        }
	    }
      }
    }
  }


  Simulation :: Simulation(Model* model, unsigned int number_trajectories, unsigned int number_timesteps, double end_time, IPropensityFunction* propensity_function, int random_seed) : model(model), end_time(end_time), random_seed(random_seed), number_timesteps(number_timesteps), number_trajectories(number_trajectories), propensity_function(propensity_function){
    timeline = new double[number_timesteps];
    double timestep_size = end_time/(number_timesteps-1);
    for(unsigned int i = 0; i < number_timesteps; i++){
      timeline[i] = timestep_size * i;
    }
    unsigned int trajectory_size = number_timesteps * (model -> number_species);
    trajectories_1D = new unsigned int[number_trajectories * trajectory_size];
    trajectories = new unsigned int**[number_trajectories];
    for(unsigned int i = 0; i < number_trajectories; i++){
      trajectories[i] = new unsigned int*[number_timesteps];
      for(unsigned int j = 0; j < number_timesteps; j++){
	trajectories[i][j] = &(trajectories_1D[i * trajectory_size + j *  (model -> number_species)]);
      }
    }    
  }


  Simulation :: ~Simulation(){
    delete timeline;
    delete trajectories_1D;
    for(unsigned int i = 0; i < number_trajectories; i++){
      delete trajectories[i];
    }
    delete trajectories;
  }

  
  std :: ostream& operator<<(std :: ostream& os, const Simulation& simulation){
    for(unsigned int i = 0; i < simulation.number_timesteps; i++){
      os << simulation.timeline[i] << " ";
      for(unsigned int trajectory = 0; trajectory < simulation.number_trajectories; trajectory++){
	for(unsigned int j = 0; j < simulation.model -> number_species; j++){
	  os << simulation.trajectories[trajectory][i][j] <<  " ";
	}
      }
      os << "\n";
    }
    return os;
  }

  void Simulation :: output_results_buffer(std :: ostream& os){
    double temp;
    unsigned char* temp_byte = reinterpret_cast<unsigned char*>(&temp);
    for(unsigned int i = 0; i < number_timesteps; i++){
      temp = timeline[i];
      for(unsigned int byte_i = 0; byte_i < sizeof(double); byte_i++){
	os << temp_byte[byte_i];
      }
      for(unsigned int trajectory = 0; trajectory < number_trajectories; trajectory++){
	for(unsigned int j = 0; j < model -> number_species; j++){
	  temp = trajectories[trajectory][i][j];
	  for(unsigned int byte_i = 0; byte_i < sizeof(double); byte_i++){
	    os << temp_byte[byte_i];
	  }
	}
      }
    }
  }
}
