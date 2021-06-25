#include <iostream>
#include "HybridModel.h"
#include "tau.h"
#include "model.h"

// toggle_reactions()
namespace Gillespy::TauHybrid {

	// Helper method to flag reactions that can be processed deterministically (continuous change)
	// without exceeding the user-supplied tolerance
	std::set<int> flag_det_rxns(
		std::vector<HybridReaction> &reactions,
		std::vector<HybridSpecies> &species)
	{
		int num_reactions = reactions.size();
		int num_species = species.size();
		std::set<int> det_rxns;

		for (int rxn_i = 0; rxn_i < reactions.size(); ++rxn_i) {
			// start with the assumption that reaction is determinstic
			HybridReaction &rxn = reactions[rxn_i];
			rxn.mode = SimulationState::CONTINUOUS;

			// iterate through the dependent species of this reaction
			// Loop breaks if we've already determined that it is to be marked as discrete.
			for (int spec_i = 0; spec_i < num_species && rxn.mode == SimulationState::CONTINUOUS; ++spec_i) {
				// Reaction has a dependency on a species if its dx is positive or negative.
				// Any species with "dependency" change of 0 is by definition not a dependency.
				if (rxn.base_reaction->species_change[spec_i] == 0) {
					continue;
				}

				// if any of the dependencies are set by the user as discrete OR
				// have been set as dynamic and has not been flagged as deterministic,
				// allow it to be modelled discretely
				if (species[spec_i].user_mode == SimulationState::DYNAMIC) {
					rxn.mode = species[spec_i].partition_mode;
				}
				else {
					rxn.mode = species[spec_i].user_mode;
				}
			}

			if (rxn.mode == SimulationState::CONTINUOUS) {
				det_rxns.insert(rxn_i);
			}
		}

		return det_rxns;
	}

	void partition_species(
		std::vector<HybridSpecies> &species,
		std::vector<HybridReaction> &reactions,
		const std::vector<double> &propensity_values, 
		std::vector<double> curr_state, 
		double tau_step, 
		const TauArgs &tauArgs)
	{
		// coefficient of variance- key:species id, value: cv
		std::map<int, double> cv;
		// means
		std::map<int, double> means;
		// standard deviation
		std::map<int, double> sd;

		// Initialize means and sd's
		for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
			HybridSpecies &spec = species[spec_i];

			if (spec.user_mode == SimulationState::DYNAMIC) {
				means.insert({ spec_i, curr_state[spec_i] });
				sd.insert({ spec_i, 0 });
			}
		}

		// calculate means and standard deviations for dynamic-mode species involved in reactions
		for (int rxn_i = 0; rxn_i < reactions.size(); ++rxn_i) {
			HybridReaction &rxn = reactions[rxn_i];

			for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
				// 0-dx species are not dependencies of this reaction, so dx == 0 is ignored.
				int spec_dx = rxn.base_reaction->species_change[spec_i];
				if (spec_dx > 0) {
					// Selected species is a reactant.
					HybridSpecies &reactant = species[spec_i];
					// Only dynamic species need to be considered.
					if (reactant.partition_mode != SimulationState::DYNAMIC) {
						continue;
					}
					means[spec_i] -= (tau_step * propensity_values[rxn_i] * spec_dx);
					sd[spec_i] += std::pow((tau_step * propensity_values[rxn_i] * spec_dx), 2);
				}
				else if (spec_dx < 0) {
					// Selected species is a product.
					HybridSpecies &product = species[spec_i];
					// Only dynamic species need to be considered.
					if (product.partition_mode != SimulationState::DYNAMIC) {
						continue;
					}
					means[spec_i] += (tau_step * propensity_values[rxn_i] * spec_dx);
					sd[spec_i] += std::pow((tau_step * propensity_values[rxn_i] * spec_dx), 2);
				}
			}
		}

		// calculate coefficient of variation using means and sd
		for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
			if (means.count(spec_i) <= 0) {
				continue;
			}

			HybridSpecies &spec = species[spec_i];
			if (spec.switch_min == 0)
			{
				// (default value means switch  min not set, use switch tol)
				if (means[spec_i] > 0) {
					cv[spec_i] = (sd[spec_i] /means[spec_i]);
				}
				else {
					cv[spec_i] = 1;
				}

				if (cv[spec_i] < spec.switch_tol) {
					spec.partition_mode = SimulationState::CONTINUOUS;
				}
				else {
					spec.partition_mode = SimulationState::DISCRETE;
				}
			}
			else
			{
				if (means[spec_i] > spec.switch_min) {
					spec.partition_mode = SimulationState::CONTINUOUS;
				}
				else {
					spec.partition_mode = SimulationState::DISCRETE;
				}
			}
		}
	}
}
