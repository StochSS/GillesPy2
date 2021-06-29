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

	HybridSimulation::HybridSimulation(const Model &model)
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

}
