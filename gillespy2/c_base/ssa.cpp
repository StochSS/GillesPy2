#include "ssa.h"
#include <random>
//Included for memcpy
#include <string.h>

namespace Gillespy{
  void ssa_direct(Simulation* simulation){
    if(simulation){
      double* propensity_values = new double[(simulation -> model) -> number_reactions];
      //Number of bytes for copying states
      int state_size = sizeof(int)*((simulation -> model) -> number_species);
      //Current state
      int* current_state = new int[(simulation -> model) -> number_species];
      //Simulate for each trajectory
      for(int trajectory_number = 0; trajectory_number < simulation -> number_trajectories; trajectory_number++){
	//Get simpler reference to memory space for this trajectory
	int** trajectory = simulation -> trajectories[trajectory_number]; 
	//Set up initial state and next_state to be written
	for(int species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
	  trajectory[0][species_number] = (simulation -> model) -> species[species_number].initial_population;
	  current_state[species_number] = trajectory[0][species_number];
	}
	double current_time = 0;
	int entry_count = 1;
	//calculate initial propensities
	for(int reaction_number = 0; reaction_number < ((simulation -> model) -> number_reactions); reaction_number++){
	  propensity_values[reaction_number] = (simulation -> propensity_function) -> evaluate(reaction_number, trajectory[0]);
	}
	double propensity_sum;
	while(current_time < (simulation -> end_time)){
	  //Sum propensities
	  propensity_sum = 0;
	  for(int reaction_number = 0; reaction_number < ((simulation -> model) -> number_reactions); reaction_number++){
	    propensity_sum += propensity_values[reaction_number];
	  }
	  //No more reactions
	  if(propensity_sum <= 0){
	    //Copy all of last changed state for rest of entries
	    for(int i = entry_count + 1; i < simulation -> number_timesteps; i++){
	      memcpy(trajectory[i], current_state, state_size);
	    }
	    //Quit simulating this trajectory
	    break;
	  }
	  //Reaction will fire, determine which one
	  double cumulative_sum = 0;//random.uniform from 0 to propensity_sum
	  current_time += 0;//-log(random.uniform)/propensity_sum
	  //Copy current state to passed timesteps
	  while(entry_count < simulation -> number_timesteps && (simulation -> timeline[entry_count]) <= current_time){
	    memcpy(trajectory[entry_count], current_state, state_size);
	    entry_count++;
	  }
	  for(int potential_reaction = 0; potential_reaction < ((simulation -> model) -> number_reactions); potential_reaction++){
	    cumulative_sum -= propensity_values[potential_reaction];
	    //This reaction fired
	    if (cumulative_sum <= 0){
	      //Update current state
	      Reaction reaction = ((simulation -> model) -> reactions[potential_reaction]);
	      for(int species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
		current_state[species_number] += reaction.species_change[species_number];
	      }
	      //Recalculate needed propensities
	      for(int& affected_reaction : reaction.affected_reactions){
		propensity_values[affected_reaction] =  (simulation -> propensity_function) -> evaluate(affected_reaction, current_state);
	      }
	      break;
	    }
	  }
	}//Simulation has reached end time
	delete propensity_values;
      }//Finished simulating all trajectories
    }//end if simulation pointer not null
  }//end ssa_direct
}//end namespace
