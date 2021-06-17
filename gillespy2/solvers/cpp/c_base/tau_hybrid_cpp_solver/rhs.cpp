#include "rhs.h"
#include "HybridModel.h"
#include "integrator.h"

using namespace Gillespy::TauHybrid;

/**
 * Integrator function for ODE linear solver.
 * This gets passed directly to the Sundials ODE solver once initialized.
 */
int Gillespy::TauHybrid::rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	// Get y(t) vector and f(t, y) vector
	realtype *Y = N_VGetArrayPointer(y);
	realtype *dydt = N_VGetArrayPointer(ydot);
	realtype propensity;

	// Extract simulation data
	IntegratorData *data = static_cast<IntegratorData*>(user_data);
	HybridSimulation *sim = data->simulation;
	HybridSpecies *species = data->species_state;
	HybridReaction *reactions = data->reaction_state;
	std::vector<double> &propensities = data->propensities;
	// Concentrations and reactions are both used for their respective propensity evaulations.
	// They both should, roughly, reflect the same data, but tau selection requires both.
	std::vector<double> &concentrations = data->concentrations;
	std::vector<int> &populations = data->populations;
	unsigned int num_species = sim->model->number_species;
	unsigned int num_reactions = sim->model->number_reactions;

	// Differentiate different regions of the input/output vectors.
	// First half is for concentrations, second half is for reaction offsets.
	realtype *rxn_offsets = &Y[num_species];
	realtype *dydt_offsets = &dydt[num_species];
	int rxn_offset_boundary = num_species + num_reactions;

	// Populate the current ODE state into the concentrations vector.
	// dy/dt results are initialized to zero, and become the change in propensity.
	unsigned int spec_i;
	for (spec_i = 0; spec_i < num_species; ++spec_i) {
		concentrations[spec_i] = Y[spec_i];
		populations[spec_i] = Y[spec_i];
		dydt[spec_i] = 0;
	}

	// Populate the current stochastic state into the root offset vector.
	// dy/dt results are initialized to zero, and become the change in offset.
	unsigned int rxn_i;
	for (rxn_i = 0; rxn_i < num_reactions; ++rxn_i) {
		dydt_offsets[rxn_i] = 0;
	}

	// Each species has a "spot" in the y and f(y,t) vector.
	// For each species, place the result of f(y,t) into dydt vector.
	int species_change;
	Gillespy::Reaction *current_rxn;

	// Process deterministic propensity state
	// These updates get written directly to the integrator's concentration state
	for (rxn_i = 0; rxn_i < num_reactions; ++rxn_i) {
		current_rxn = reactions[rxn_i].base_reaction;
		// NOTE: we may need to evaluate ODE and Tau propensities separately.
		// At the moment, it's unsure whether or not that's required.

		switch (reactions[rxn_i].mode) {
		case SimulationState::DISCRETE:
			// Process stochastic reaction state by updating the root offset for each reaction.
			propensity = sim->propensity_function->TauEvaluate(rxn_i, populations);
			dydt_offsets[rxn_i] += propensity;
			// break; left out on purpose, continuous results happen no matter what
		case SimulationState::CONTINUOUS:
			propensity = sim->propensity_function->ODEEvaluate(rxn_i, concentrations);
			propensities[rxn_i] = propensity;

			for (spec_i = 0; spec_i < num_species; ++spec_i) {
				// Use the evaluated propensity to update the concentration levels and reaction state.
				// Propensity is treated as positive if it's a product, negative if it's a reactant.
				species_change = current_rxn->species_change[spec_i];
				if (species_change == 0)
					continue;

				// The product on the right evaluates to 1 if species_change is positive,
				//    and -1 if it's negative.
				// This is a branchless alternative to using an if-statement.
				dydt[spec_i] += propensity * (-1 + 2 * (species_change > 0));
			}
			break;
		default:
			break;
		}
	}

	return 0;
};
