#pragma once

#include <variant>
#include "model.h"

#define GPY_HYBRID_ABSTOL 1e-5
#define GPY_HYBRID_RELTOL 1e-5

namespace Gillespy::TauHybrid {

	enum SimulationState
	{
		CONTINUOUS = 0,
		DISCRETE,
		DYNAMIC
	};

	struct HybridSpecies
	{
		Species *base_species;

		// allows the user to specify if a species' population should definitely be modeled continuously or 
		// discretely
		// CONTINUOUS or DISCRETE
		// otherwise, mode will be determined by the program (DYNAMIC)
		// if no choice is made, DYNAMIC will be assumed 
		SimulationState  user_mode : 2;

		// during simulation execution, a species will fall into either of the two categories, CONTINUOUS or DISCRETE
		// this is pre-determined only if the user_mode specifies CONTINUOUS or DISCRETE.
		// otherwise, if DYNAMIC is specified, partition_mode will be continually calculated throughout the simulation
		// according to standard deviation and coefficient of variance.
		SimulationState partition_mode : 1;

		// Tolerance level for considering a dynamic species deterministically, value is compared
		// to an estimated sd/mean population of a species after a given time step.
		//  This value will be used if a switch_min is not provided. The default value is 0.03
		double switch_tol = 0.03;

		//Minimum population value at which species will be represented as continuous. 
		// If a value is given, switch_min will be used instead of switch_tol.
		unsigned int switch_min = 0;

		HybridSpecies();
	};

	struct HybridReaction
	{
		Reaction *base_reaction;
		SimulationState mode : 1;

		HybridReaction();
	};

	union hybrid_state
	{
		unsigned int discrete;
		double continuous;
	};

	struct HybridSimulation : Simulation<double>
	{
		hybrid_state *trajectories_hybrid1D;
		hybrid_state ***trajectories_hybrid;

		double hybrid_propensity(int reaction_id, std::vector<double> &current_state);

		HybridSimulation();
	};

}
