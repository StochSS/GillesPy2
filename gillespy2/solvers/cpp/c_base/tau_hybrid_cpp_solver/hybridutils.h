#pragma once

#include "HybridModel.h"
#include "tau.h"

namespace Gillespy::TauHybrid
{

	std::set<int> flag_det_rxns(
		std::vector<HybridReaction> &reactions,
		std::vector<HybridSpecies> &species);

	void partition_species(
		std::vector<HybridReaction> &reactions,
		std::vector<HybridSpecies> &species,
		const std::vector<double> &propensity_values,
		std::vector<double> &curr_state,
		double tau_step,
		const TauArgs &TauArgs);

}
