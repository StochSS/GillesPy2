#include "ssa.h"
#include <cmath>//Included for natural logarithm
#include <string.h>//Included for memcpy only

namespace Gillespy{
  void ode(Simulation* simulation){
    if(simulation){
      //Number of bytes for copying states
      unsigned int state_size = sizeof(int)*((simulation -> model) -> number_species);
      //Current state
      unsigned int* current_state = new unsigned int[(simulation -> model) -> number_species];
      //Calculated propensity values for current state
      double* propensity_values = new double[(simulation -> model) -> number_reactions];

      //copy initial state for each trajectory
      for(unsigned int species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
	simulation -> trajectories[0][0][species_number] = (simulation -> model) -> species[species_number].initial_population;
      }
      //Simulate for each trajectory
      for(unsigned int trajectory_number = 0; trajectory_number < simulation -> number_trajectories; trajectory_number++){
	//Get simpler reference to memory space for this trajectory
	unsigned int** trajectory = simulation -> trajectories[trajectory_number];
	//Copy initial state as needed
	if(trajectory_number > 0){
	  memcpy(trajectory[0], simulation -> trajectories[0][0], state_size);
	}
	//Set up current state from initial state
	memcpy(current_state, trajectory[0], state_size);
	double current_time = 0;
	unsigned int entry_count = 1;
	//calculate initial propensities
	for(unsigned int reaction_number = 0; reaction_number < ((simulation -> model) -> number_reactions); reaction_number++){
	  propensity_values[reaction_number] = (simulation -> propensity_function) -> evaluate(reaction_number, current_state);
	}
	double propensity_sum;
	while(current_time < (simulation -> end_time)){

	}//Simulation has reached end time
      }//Finished simulating all trajectories
      delete propensity_values;
      delete current_state;
    }//end if simulation pointer not null
  }//end ode
}//end namespace
