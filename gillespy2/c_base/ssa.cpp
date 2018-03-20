#include "ssa.h"
#include <random>//Included for mt19937 random number generator
#include <cmath>//Included for natural logarithm
#include <string.h>//Included for memcpy only
#include <iostream>

namespace Gillespy{
  void ssa_direct(Simulation* simulation){
    if(simulation){
      std :: mt19937_64 rng(simulation -> random_seed);
      //Number of bytes for copying states
      uint state_size = sizeof(int)*((simulation -> model) -> number_species);
      //Current state
      uint* current_state = new uint[(simulation -> model) -> number_species];
      //Calculated propensity values for current state
      double* propensity_values = new double[(simulation -> model) -> number_reactions];
      
      //copy initial state for each trajectory
      for(uint species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
	simulation -> trajectories[0][0][species_number] = (simulation -> model) -> species[species_number].initial_population;
      }
      //Simulate for each trajectory
      for(uint trajectory_number = 0; trajectory_number < simulation -> number_trajectories; trajectory_number++){
	//std :: cout << "Simulating trajectory " << trajectory_number << std :: endl;
	//Get simpler reference to memory space for this trajectory
	uint** trajectory = simulation -> trajectories[trajectory_number];
	//Copy initial state as needed
	if(trajectory_number > 0){
	  memcpy(trajectory[0], simulation -> trajectories[0][0], state_size);
	}
	memcpy(current_state, trajectory[0], state_size);
	//Set up initial state and next_state to be written
	for(uint species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
	  trajectory[0][species_number] = (simulation -> model) -> species[species_number].initial_population;
	}
	double current_time = 0;
	uint entry_count = 1;
	//calculate initial propensities
	for(uint reaction_number = 0; reaction_number < ((simulation -> model) -> number_reactions); reaction_number++){
	  propensity_values[reaction_number] = (simulation -> propensity_function) -> evaluate(reaction_number, trajectory[0]);
	}
	double propensity_sum;
	while(current_time < (simulation -> end_time)){
	  /* std :: cout << "Current time: " << current_time << "\t Entries: " << entry_count << ": ";
	  for(uint i = 0; i < simulation -> model -> number_species; i++){
	    std :: cout << current_state[i] << ", ";
	  }
	  std :: cout << std :: endl;*/
	  //Sum propensities
	  propensity_sum = 0;
	  for(uint reaction_number = 0; reaction_number < ((simulation -> model) -> number_reactions); reaction_number++){
	    propensity_sum += propensity_values[reaction_number];
	  }
	  //No more reactions
	  if(propensity_sum <= 0){
	    //Copy all of last changed state for rest of entries
	    for(uint i = entry_count + 1; i < simulation -> number_timesteps; i++){
	      memcpy(trajectory[i], current_state, state_size);
	    }
	    //Quit simulating this trajectory
	    break;
	  }//End if no more reactions
	  
	  //Reaction will fire, determine which one
	  double cumulative_sum = rng() * propensity_sum/rng.max();//random.uniform from 0 to propensity_sum
	  current_time += -log(rng() * 1.0 / rng.max()) / propensity_sum;//-log(random.uniform)/propensity_sum
	  //Copy current state to passed timesteps
	  while(entry_count < simulation -> number_timesteps && (simulation -> timeline[entry_count]) <= current_time){
	    memcpy(trajectory[entry_count], current_state, state_size);
	    entry_count++;
	  }
	  
	  for(uint potential_reaction = 0; potential_reaction < ((simulation -> model) -> number_reactions); potential_reaction++){
	    cumulative_sum -= propensity_values[potential_reaction];
	    //This reaction fired
	    if (cumulative_sum <= 0 && propensity_values[potential_reaction] > 0){
	      //Update current state
	      Reaction reaction = ((simulation -> model) -> reactions[potential_reaction]);
	      for(uint species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
		current_state[species_number] += reaction.species_change[species_number];
	      }
	      //Recalculate needed propensities
	      for(uint& affected_reaction : reaction.affected_reactions){
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
