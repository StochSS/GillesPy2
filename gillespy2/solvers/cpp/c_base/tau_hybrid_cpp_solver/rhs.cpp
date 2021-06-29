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
	std::vector<HybridSpecies> *species = data->species_state;
	std::vector<HybridReaction> *reactions = data->reaction_state;
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
	}

	// Deterministic reactions generally are "evaluated" by generating dy/dt functions
	//   for each of their dependent species.
	// To handle these, we will go ahead and evaluate each species' differential equations.
	for (spec_i = 0; spec_i < num_species; ++spec_i) {
		dydt[spec_i] = (*species)[spec_i].diff_equation.evaluate(concentrations, populations);
	}

	// Process deterministic propensity state
	// These updates get written directly to the integrator's concentration state
	for (int rxn_i = 0; rxn_i < num_reactions; ++rxn_i) {
		switch ((*reactions)[rxn_i].mode) {
		case SimulationState::DISCRETE:
			// Process stochastic reaction state by updating the root offset for each reaction.
			propensity = (*reactions)[rxn_i].ssa_propensity(rxn_i, populations);
			dydt_offsets[rxn_i] = propensity;
			propensities[rxn_i] = propensity;
			break;

		case SimulationState::CONTINUOUS:
		default:
			dydt_offsets[rxn_i] = 0;
			break;
		}
	}

	return 0;
};
