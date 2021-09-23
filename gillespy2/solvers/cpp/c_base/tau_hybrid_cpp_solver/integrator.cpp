/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2021 GillesPy2 developers.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "integrator.h"

using namespace Gillespy::TauHybrid;

static bool validate(Integrator *integrator, int retcode);

IntegratorData::IntegratorData(
	HybridSimulation *simulation,
	int num_species,
	int num_reactions)
	: simulation(simulation),
	  concentrations(std::vector<double>(num_species)),
	  populations(std::vector<int>(num_species)),
	  propensities(std::vector<double>(num_reactions)),
	  species_state(&simulation->species_state),
	  reaction_state(&simulation->reaction_state) {}

IntegratorData::IntegratorData(HybridSimulation *simulation)
	: IntegratorData(
		simulation,
		simulation->model->number_species,
		simulation->model->number_reactions) {}


Integrator::Integrator(HybridSimulation *simulation, N_Vector y0, double reltol, double abstol)
	: y0(y0),
	  t(0.0f),
	  t0(0.0f),
	  y(N_VClone_Serial(y0)),
	  data(simulation),
	  num_reactions(simulation->model->number_reactions),
	  num_species(simulation->model->number_species)
{
	// y0 is the initial state, y is updated during integration.
	// N_VClone_Serial() does not clone *contents*, we have to do that explicitly.
	for (int mem_i = 0; mem_i < num_reactions + num_species; ++mem_i) {
		NV_Ith_S(y, mem_i) = NV_Ith_S(this->y0, mem_i);
	}

	for (int spec_i = 0; spec_i < num_species; ++spec_i)
	{
		data.populations[spec_i]
			= data.concentrations[spec_i]
			= simulation->model->species[spec_i].initial_population;
	}

	for (int rxn_i = 0; rxn_i < num_reactions; ++rxn_i)
	{
		data.propensities[rxn_i] = 0;
	}

	cvode_mem = CVodeCreate(CV_BDF);
	validate(this, CVodeInit(cvode_mem, rhs, t, y));
	validate(this, CVodeSStolerances(cvode_mem, reltol, abstol));

	solver = SUNLinSol_SPGMR(y, 0, 0);
	validate(this, CVodeSetUserData(cvode_mem, &data));
	validate(this, CVodeSetLinearSolver(cvode_mem, solver, NULL));
}

double Integrator::save_state()
{
	int max_offset = num_reactions + num_species;
	for (int mem_i = 0; mem_i < max_offset; ++mem_i)
	{
		NV_Ith_S(y0, mem_i) = NV_Ith_S(y, mem_i);
	}

	t0 = t;
	return t0;
}

double Integrator::restore_state()
{
	int max_offset = num_reactions + num_species;
	for (int mem_i = 0; mem_i < max_offset; ++mem_i)
	{
		NV_Ith_S(y, mem_i) = NV_Ith_S(y0, mem_i);
	}
	if (!validate(this, CVodeReInit(cvode_mem, t0, y0)))
	{
		return 0;
	}

	return t0;
}

void Integrator::refresh_state()
{
	validate(this, CVodeReInit(cvode_mem, t, y));
}

void Integrator::reinitialize(N_Vector y_reset)
{
	int max_offset = num_reactions + num_species;
	for (int mem_i = 0; mem_i < max_offset; ++mem_i)
	{
		NV_Ith_S(y0, mem_i) = NV_Ith_S(y_reset, mem_i);
	}
	validate(this, CVodeReInit(cvode_mem, 0, y0));
}

Integrator::~Integrator()
{
	N_VDestroy_Serial(y);
	CVodeFree(&cvode_mem);
	SUNLinSolFree_SPGMR(solver);
}

IntegrationResults Integrator::integrate(double *t)
{
	if (!validate(this, CVode(cvode_mem, *t, y, &this->t, CV_NORMAL)))
	{
		return { nullptr, nullptr };
	}

	return {
		NV_DATA_S(y), // NV_DATA_S instead?
		NV_DATA_S(y) + num_species
	};
}


URNGenerator::URNGenerator()
	: uniform(0, 1) {}

URNGenerator::URNGenerator(unsigned long long seed)
	: uniform(0, 1),
	  rng(seed)
{
	this->seed = seed;
}


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
N_Vector Gillespy::TauHybrid::init_model_vector(Gillespy::Model<double> &model, URNGenerator urn)
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
	for (int spec_i = 0; spec_i < model.number_species; ++spec_i)
	{
		NV_Ith_S(y0, spec_i) = model.species[spec_i].initial_population;
	}

	// The second half represents the current "randomized state" for each reaction.
	// ... | --- rxn_offsets --- ]
	for (int rxn_i = model.number_species; rxn_i < rxn_offset_boundary; ++rxn_i)
	{
		// Represents the current "randomized state" for each reaction, used as a
		//   helper value to determine if/how many stochastic reactions fire.
		// This gets initialized to a random negative offset, and gets "less negative"
		//   during the integration step.
		// After each integration step, the reaction_state is used to count stochastic reactions.
		NV_Ith_S(y0, rxn_i) = log(urn.next());
	}

	return y0;
}

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
	std::vector<int> &populations = data->populations;
	unsigned int num_species = sim->model->number_species;
	unsigned int num_reactions = sim->model->number_reactions;

	// Differentiate different regions of the input/output vectors.
	// First half is for concentrations, second half is for reaction offsets.
	realtype *dydt_offsets = &dydt[num_species];

	// Populate the current ODE state into the concentrations vector.
	// dy/dt results are initialized to zero, and become the change in propensity.
	unsigned int spec_i;
	for (spec_i = 0; spec_i < num_species; ++spec_i)
	{
		populations[spec_i] = static_cast<int>(Y[spec_i]);
	}

	// Deterministic reactions generally are "evaluated" by generating dy/dt functions
	//   for each of their dependent species.
	// To handle these, we will go ahead and evaluate each species' differential equations.
	for (spec_i = 0; spec_i < num_species; ++spec_i)
	{
		if ((*species)[spec_i].boundary_condition) {
			// The effective dy/dt of a boundary condition is 0.
			dydt[spec_i] = 0.0;
		}
		else
		{
			dydt[spec_i] = (*species)[spec_i].diff_equation.evaluate(t, Y, &populations[0]);
		}
	}

	// Process deterministic propensity state
	// These updates get written directly to the integrator's concentration state
	for (int rxn_i = 0; rxn_i < num_reactions; ++rxn_i)
	{
		switch ((*reactions)[rxn_i].mode) {
		case SimulationState::DISCRETE:
			// Process stochastic reaction state by updating the root offset for each reaction.
			propensity = HybridReaction::ssa_propensity(rxn_i, &populations[0]);
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


static bool validate(Integrator *integrator, int retcode)
{
	switch (retcode)
	{
	case CV_MEM_NULL:
		integrator->status = IntegrationStatus::NULL_POINTER;
		return false;
	case CV_NO_MALLOC:
		integrator->status = IntegrationStatus::BAD_MEMORY;
		return false;
	case CV_TOO_CLOSE:
	case CV_TOO_MUCH_WORK:
		integrator->status = IntegrationStatus::BAD_STEP_SIZE;
		return false;
	case CV_SUCCESS:
	default:
		integrator->status = IntegrationStatus::OK;
		return true;
	}
}
