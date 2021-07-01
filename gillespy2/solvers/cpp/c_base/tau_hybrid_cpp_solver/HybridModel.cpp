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
		: user_mode(SimulationState::CONTINUOUS),
		  partition_mode(SimulationState::CONTINUOUS),
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

	HybridSimulation::HybridSimulation(const Model<double> &model)
	    : Simulation<double>(),
		  species_state(model.number_species),
		  reaction_state(model.number_reactions)
    {
		for (int spec_i = 0; spec_i < model.number_species; ++spec_i)
		{
			species_state[spec_i].base_species = &model.species[spec_i];
		}

		for (int rxn_i = 0; rxn_i < model.number_reactions; ++rxn_i)
		{
			reaction_state[rxn_i].base_reaction = &model.reactions[rxn_i];
		}
    }


	double DifferentialEquation::evaluate(
		std::vector<double> &ode_state,
		std::vector<int> &ssa_state)
	{
		double sum =  0.0;

		for (auto &formula : formulas) {
			sum += formula(ode_state, ssa_state);
		}

		return sum;
	}


	void create_differential_equations(
		std::vector<HybridSpecies> &species,
		std::vector<HybridReaction> &reactions)
	{
		// For now, differential equations are generated from scratch.
		// It may be more efficient to determine which formulas need to change.
		// Until then, the compound formulas in every species are cleared.
		for (HybridSpecies &spec : species) {
			spec.diff_equation.formulas.clear();
		}

		for (int rxn_i = 0; rxn_i < reactions.size(); ++rxn_i) {
			HybridReaction rxn = reactions[rxn_i];
			if (rxn.mode == SimulationState::DISCRETE) {
				continue;
			}

			for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
				// A species change of 0 indicates that this species is not a dependency for this reaction.
				if (rxn.base_reaction->species_change[spec_i] == 0) {
					continue;
				}

				HybridSpecies &spec = species[spec_i];
				auto &formula_set = spec.diff_equation.formulas;
				int spec_diff = rxn.base_reaction->species_change[spec_i];

				switch (spec.partition_mode) {
				case SimulationState::CONTINUOUS:
					formula_set.push_back([&rxn, rxn_i, spec_diff](
						std::vector<double> &ode_state,
						std::vector<int> &ssa_state)
					{
						return spec_diff * rxn.ode_propensity(rxn_i, ode_state);
					});
					break;

				case SimulationState::DISCRETE:
					spec.diff_equation.formulas.push_back([&rxn, rxn_i, spec_diff](
						std::vector<double> &ode_state,
						std::vector<int> &ssa_state)
					{
						return spec_diff * rxn.ssa_propensity(rxn_i, ssa_state);
					});
					break;

				default:
					break;
				}
			}
		}
	}

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
		std::vector<HybridReaction> &reactions,
		std::vector<HybridSpecies> &species,
		const std::vector<double> &propensity_values, 
		std::vector<double> &curr_state, 
		double tau_step, 
		const TauArgs<double> &tauArgs)
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
				// Only dynamic species whose mean/SD is requested are to be considered.
				if (means.count(spec_i) <= 0) {
					continue;
				}
				// Selected species is either a reactant or a product, depending on whether
				//   dx is positive or negative.
				// 0-dx species are not dependencies of this reaction, so dx == 0 is ignored.
				int spec_dx = rxn.base_reaction->species_change[spec_i];
				if (spec_dx < 0) {
					// Selected species is a reactant.
					means[spec_i] -= (tau_step * propensity_values[rxn_i] * spec_dx);
					sd[spec_i] += std::pow((tau_step * propensity_values[rxn_i] * spec_dx), 2);
				}
				else if (spec_dx > 0) {
					// Selected species is a product.
					HybridSpecies &product = species[spec_i];
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

				spec.partition_mode = cv[spec_i] < spec.switch_tol
					? SimulationState::CONTINUOUS
					: SimulationState::DISCRETE;
			}
			else
			{
				spec.partition_mode = means[spec_i] > spec.switch_min
					? SimulationState::CONTINUOUS
					: SimulationState::DISCRETE;
			}
		}
	}

	void update_species_state(
		std::vector<HybridSpecies> &species,
		std::vector<double> &current_state)
	{
		for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
			switch (species[spec_i].partition_mode) {
			case SimulationState::CONTINUOUS:
				current_state[spec_i] = std::floor(current_state[spec_i]);
				break;
			default:
				break;
			}
		}
	}

}
