#include "HybridModel.h"

namespace Gillespy::TauHybrid
{

	HybridReaction::HybridReaction()
		: mode(SimulationState::DISCRETE),
		  base_reaction(nullptr)
	{
		// Empty constructor body
	}

	HybridSpecies::HybridSpecies()
		: user_mode(SimulationState::DISCRETE),
		  partition_mode(SimulationState::DISCRETE),
		  switch_tol(0.03),
		  switch_min(0)
	{
		// Empty constructor body
	}

	HybridSimulation::HybridSimulation()
		: Simulation<double>()
	{
		// Empty constructor body
	}

}
