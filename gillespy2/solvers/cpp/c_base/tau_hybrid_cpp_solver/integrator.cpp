#include "integrator.h"
#include "rhs.h"

static bool validate(int retcode);

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
		simulation->model->number_reactions) {}

IntegratorData::~IntegratorData()
{
	delete[] species_state;
	delete[] reaction_state;
}


Integrator::Integrator(HybridSimulation *simulation, N_Vector y0, double reltol, double abstol)
	: y(y0),
	  data(simulation),
	  num_reactions(simulation->model->number_reactions),
	  num_species(simulation->model->number_species)
{
	for (int spec_i = 0; spec_i < num_species; ++spec_i) {
		data.populations[spec_i]
			= data.concentrations[spec_i]
			= simulation->model->species[spec_i].initial_population;
	}

	cvode_mem = CVodeCreate(CV_BDF);
	validate(CVodeInit(cvode_mem, rhs, t, y));
	validate(CVodeSStolerances(cvode_mem, reltol, abstol));

	solver = SUNLinSol_SPGMR(y, 0, 0);
	validate(CVodeSetUserData(cvode_mem, &data));
	validate(CVodeSetLinearSolver(cvode_mem, solver, NULL));
}

Integrator::~Integrator()
{
	N_VDestroy_Serial(y);
	CVodeFree(&cvode_mem);
	SUNLinSolFree_SPGMR(solver);
}

IntegrationResults Integrator::integrate(double *t)
{
	if (!validate(CVode(cvode_mem, *t, y, t, CV_NORMAL))) {
		return { nullptr, nullptr };
	}

	return {
		NV_DATA_S(y), // NV_DATA_S instead?
		NV_DATA_S(y) + num_species
	};
}


URNGenerator::URNGenerator()
	: uniform(0, 1) {}

URNGenerator::URNGenerator(double seed)
	: uniform(0, 1),
	  rng(seed) {}


/* Generate a new random floating-point number on the range [0,1).
 * Uses a uniform distribution to generate.
 */
double URNGenerator::next()
{
	return uniform(rng);
}


/* Initialize a SUNDials N_Vector based on information provided in the model.
 * 
 */
N_Vector Gillespy::TauHybrid::init_model_vector(Gillespy::Model &model, URNGenerator urn)
{
	int rxn_offset_boundary = model.number_reactions + model.number_species;

	// INITIAL INTEGRATOR STATE VECTOR
	// Integrator is used to integrate two vector regions separately:
	//   - concentrations for deterministic reactions
	//   - reaction offsets for stochastic reactions
	// [ --- concentrations --- | --- rxn_offsets --- ]
	// concentrations: bounded by [0, num_species)
	// rxn_offsets:    bounded by [num_species, num_species + num_reactions)
	N_Vector y0 = N_VNew_Serial(rxn_offset_boundary);

	// The first half of the integration vector is used for integrating species concentrations.
	// [ --- concentrations --- | ...
	for (int spec_i = 0; spec_i < model.number_species; ++spec_i) {
		NV_Ith_S(y0, spec_i) = model.species[spec_i].initial_population;
	}

	// The second half represents the current "randomized state" for each reaction.
	// ... | --- rxn_offsets --- ]
	for (int rxn_i = model.number_species; rxn_i < rxn_offset_boundary; ++rxn_i) {
		// Represents the current "randomized state" for each reaction, used as a
		//   helper value to determine if/how many stochastic reactions fire.
		// This gets initialized to a random negative offset, and gets "less negative"
		//   during the integration step.
		// After each integration step, the reaction_state is used to count stochastic reactions.
		NV_Ith_S(y0, rxn_i) = log(urn.next());
	}

	return y0;
}


bool validate(int retcode)
{
	return true;
}
