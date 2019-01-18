#include "ssa.h"
#include <random>//Included for mt19937 random number generator
#include <cmath>//Included for natural logarithm
#include <string.h>//Included for memcpy only
#include <cstdint>

namespace Gillespy{
  void ssa_direct(Simulation* simulation){
    if(simulation){
      std :: mt19937_64 rng(simulation -> random_seed);
      //Number of bytes for copying states
      uint32_t state_size = sizeof(int)*((simulation -> model) -> number_species);
      //Current state
      uint32_t* current_state = new uint32_t[(simulation -> model) -> number_species];
      //Calculated propensity values for current state
      double* propensity_values = new double[(simulation -> model) -> number_reactions];
      
      //copy initial state for each trajectory
      for(uint32_t species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
	simulation -> trajectories[0][0][species_number] = (simulation -> model) -> species[species_number].initial_population;
      }
      //Simulate for each trajectory
      for(uint32_t trajectory_number = 0; trajectory_number < simulation -> number_trajectories; trajectory_number++){
	//Get simpler reference to memory space for this trajectory
	uint32_t** trajectory = simulation -> trajectories[trajectory_number];
	//Copy initial state as needed
	if(trajectory_number > 0){
	  memcpy(trajectory[0], simulation -> trajectories[0][0], state_size);
	}
	//Set up current state from initial state
	memcpy(current_state, trajectory[0], state_size);
	double current_time = 0;
	uint32_t entry_count = 1;
	//calculate initial propensities
	for(uint32_t reaction_number = 0; reaction_number < ((simulation -> model) -> number_reactions); reaction_number++){
	  propensity_values[reaction_number] = (simulation -> propensity_function) -> evaluate(reaction_number, current_state);
	}
	double propensity_sum;
	while(current_time < (simulation -> end_time)){
	  //Sum propensities
	  propensity_sum = 0;
	  for(uint32_t reaction_number = 0; reaction_number < ((simulation -> model) -> number_reactions); reaction_number++){
	    propensity_sum += propensity_values[reaction_number];
	  }
	  //No more reactions
	  if(propensity_sum <= 0){
	    //Copy all of last changed state for rest of entries
	    for(uint32_t i = entry_count; i < simulation -> number_timesteps; i++){
	      memcpy(trajectory[i], current_state, state_size);
	    }
	    //Quit simulating this trajectory
	    break;
	  }//End if no more reactions
	  
	  //Reaction will fire, determine which one
	  double cumulative_sum = rng() * propensity_sum/rng.max();
	  current_time += -log(rng() * 1.0 / rng.max()) / propensity_sum;
	  //Copy current state to passed timesteps
	  while(entry_count < simulation -> number_timesteps && (simulation -> timeline[entry_count]) <= current_time){
	    memcpy(trajectory[entry_count], current_state, state_size);
	    entry_count++;
	  }
	  
	  for(uint32_t potential_reaction = 0; potential_reaction < ((simulation -> model) -> number_reactions); potential_reaction++){
	    cumulative_sum -= propensity_values[potential_reaction];
	    //This reaction fired
	    if (cumulative_sum <= 0 && propensity_values[potential_reaction] > 0){
	      //Update current state
	      Reaction& reaction = ((simulation -> model) -> reactions[potential_reaction]);
	      for(uint32_t species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
		current_state[species_number] += reaction.species_change[species_number];
	      }
	      //Recalculate needed propensities
	      for(uint32_t& affected_reaction : reaction.affected_reactions){
	 	propensity_values[affected_reaction] =  (simulation -> propensity_function) -> evaluate(affected_reaction, current_state);
	      }
	      break;
	    }//Finished updating state/propensities with this reaction
	  }//Finished checking for which reaction fired at this time
	}//Simulation has reached end time
      }//Finished simulating all trajectories
      delete propensity_values;
      delete current_state;
    }//end if simulation pointer not null
  }//end ssa_direct
}//end namespace
