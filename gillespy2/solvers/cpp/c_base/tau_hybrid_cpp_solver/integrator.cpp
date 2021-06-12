#include "integrator.h"

using namespace Gillespy::TauHybrid;

IntegratorData::IntegratorData(
	HybridSimulation *simulation,
	int num_species,
	int num_reactions)
	: simulation(simulation),
	  concentrations(std::vector<double>(num_species)),
	  populations(std::vector<int>(num_species)),
	  propensities(std::vector<double>(num_reactions))
{
	species_state = new HybridSpecies[num_species];
	for (int spec_i = 0; spec_i < num_species; ++spec_i) {
		species_state[spec_i].base_species = &simulation->model->species[spec_i];
	}

	reaction_state = new HybridReaction[num_reactions];
	for (int rxn_i = 0; rxn_i < num_reactions; ++rxn_i) {
		reaction_state[rxn_i].base_reaction = &simulation->model->reactions[rxn_i];
	}
}

IntegratorData::IntegratorData(HybridSimulation *simulation)
	: IntegratorData(
		simulation,
		simulation->model->number_species,
		simulation->model->number_reactions)
{
	// Empty constructor body
}

IntegratorData::~IntegratorData()
{
	delete[] species_state;
	delete[] reaction_state;
}
