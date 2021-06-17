#include "HybridModel.h"

namespace Gillespy::TauHybrid
{

	HybridReaction::HybridReaction()
		: mode(SimulationState::CONTINUOUS),
		  base_reaction(nullptr)
	{
		// Empty constructor body
	}

	HybridSpecies::HybridSpecies()
		: user_mode(SimulationState::CONTINUOUS),
		  partition_mode(SimulationState::CONTINUOUS),
		  switch_tol(0.03),
		  switch_min(0)
	{
		// Empty constructor body
	}

	HybridSimulation::~HybridSimulation()
	{
		if (type == HYBRID) {
			for(unsigned int i = 0; i < number_trajectories; i++){
				delete trajectoriesODE[i];
			}
			delete trajectoriesODE;
		}
	}

}
