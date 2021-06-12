#pragma once

#include <vector>
#include "HybridModel.h"

namespace Gillespy::TauHybrid
{

	struct IntegratorData
	{
		HybridSimulation *simulation;
		HybridSpecies *species_state;
		HybridReaction *reaction_state;

		std::vector<double> concentrations;
		std::vector<int> populations;
		std::vector<double> propensities;

		IntegratorData(HybridSimulation *simulation);
		IntegratorData(HybridSimulation *simulation, int num_species, int num_reactions);
		~IntegratorData();
	};

}
