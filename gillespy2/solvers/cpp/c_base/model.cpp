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

  void simulationSSAINIT(Model* model, Simulation &simulation){
  	simulation.timeline = new double[simulation.number_timesteps];
    double timestep_size = simulation.end_time/(simulation.number_timesteps-1);
    for(unsigned int i = 0; i < simulation.number_timesteps; i++){
    	simulation.timeline[i] = timestep_size * i;
    }
 	unsigned int trajectory_size = simulation.number_timesteps * (model -> number_species);
    simulation.trajectories_1D = new unsigned int[simulation.number_trajectories * trajectory_size];
    simulation.trajectories = new unsigned int**[simulation.number_trajectories];
    for(unsigned int i = 0; i < simulation.number_trajectories; i++){
    	simulation.trajectories[i] = new unsigned int*[simulation.number_timesteps];
    	for(unsigned int j = 0; j < simulation.number_timesteps; j++){
    		simulation.trajectories[i][j] = &(simulation.trajectories_1D[i * trajectory_size + j *  (model -> number_species)]);
      }
    }
  }


void simulationODEINIT(Model* model, Simulation &simulation){
    simulation.timeline = new double[simulation.number_timesteps];
    double timestep_size = simulation.end_time/(simulation.number_timesteps-1);
    for(unsigned int i = 0; i < simulation.number_timesteps; i++){
      simulation.timeline[i] = timestep_size * i;
    }


    unsigned int trajectory_size = simulation.number_timesteps * (model -> number_species);
    simulation.trajectories_1DODE = new double[simulation.number_trajectories * trajectory_size];
    simulation.trajectoriesODE = new double**[simulation.number_trajectories];
    for(unsigned int i = 0; i < simulation.number_trajectories; i++){
      simulation.trajectoriesODE[i] = new double*[simulation.number_timesteps];
      for(unsigned int j = 0; j < simulation.number_timesteps; j++){
	simulation.trajectoriesODE[i][j] = &(simulation.trajectories_1DODE[i * trajectory_size + j *  (model -> number_species)]);
      }
    }
  }


  Simulation :: ~Simulation(){
    delete timeline;
    if (ISODE==1){
    delete trajectories_1DODE;
    for(unsigned int i = 0; i < number_trajectories; i++){
      delete trajectoriesODE[i];
    }
    delete trajectoriesODE;
    }else{

    delete trajectories_1D;
     for(unsigned int i = 0; i < number_trajectories; i++){
      delete trajectories[i];
    }
    delete trajectories;
    }
  }


  std :: ostream& operator<<(std :: ostream& os, const Simulation& simulation){
    for(unsigned int i = 0; i < simulation.number_timesteps; i++){
      os << simulation.timeline[i] << " ";
      for(unsigned int trajectory = 0; trajectory < simulation.number_trajectories; trajectory++){
	for(unsigned int j = 0; j < simulation.model -> number_species; j++){
	    if (simulation.ISODE==1){os << simulation.trajectoriesODE[trajectory][i][j] <<  " ";}
	    else{os << simulation.trajectories[trajectory][i][j] <<  " ";}
	}
      }
      os << "\n";
    }
    return os;
  }

void Simulation :: output_results_buffer(std::ostream& os){
    for (int i = 0 ; i < number_trajectories; i++){
        for (int j = 0; j<number_timesteps;j++){
            os<<timeline[j]<<',';
            for (int k = 0; k<model->number_species; k++){
                if (ISODE==1){os<<trajectoriesODE[i][j][k]<<',';}
                else{os<<trajectories[i][j][k]<<',';}
                }
            }
         }
    os<<(int)current_time;
    }


}