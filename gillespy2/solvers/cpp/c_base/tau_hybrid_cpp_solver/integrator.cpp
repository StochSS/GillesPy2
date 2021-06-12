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

IntegrationResults Integrator::integrate(double &t)
{
	if (!validate(CVode(cvode_mem, t, y, &t, CV_NORMAL))) {
		return { nullptr, nullptr };
	}

	return {
		NV_DATA_S(y), // NV_DATA_S instead?
		NV_DATA_S(y) + num_species
	};
}

bool validate(int retcode)
{
	return true;
}
