#pragma once

#include "HybridModel.h"
#include "cvode.h"
#include "sunlinsol_spgmr.h"
#include "sundials_types.h"
#include "nvector_serial.h"
#include <vector>
#include <random>

namespace Gillespy::TauHybrid
{

	struct IntegratorData
	{
		HybridSimulation *simulation;
		HybridSpecies *species_state;
		HybridReaction *reaction_state;

		std::vector<double> concentrations;
		std::vector<int> populations;
		std::vector<double> propensities;

		IntegratorData(HybridSimulation *simulation);
		IntegratorData(HybridSimulation *simulation, int num_species, int num_reactions);
		IntegratorData(IntegratorData &prev_data);
		~IntegratorData();
	};

	// [ --- concentrations --- | --- rxn_offsets --- ]
	struct IntegrationResults
	{
		// concentrations: bounded by [0, num_species)
		realtype *concentrations;
		// reactions:      bounded by [num_species, num_species + num_reactions)
		realtype *reactions;
	};

	class Integrator
	{
	private:
		void *cvode_mem;
		SUNLinearSolver solver;
		int num_species;
		int num_reactions;
	public:
		N_Vector y;
		realtype t;
		IntegrationResults integrate(double &t);
		IntegratorData data;

		Integrator(HybridSimulation *simulation);
		Integrator(HybridSimulation *simulation, N_Vector y0, double reltol, double abstol);
		~Integrator();
	};

	struct URNGenerator
	{
	private:
		std::uniform_real_distribution<double> uniform;
		std::mt19937_64 rng;
	public:
		double next();
		URNGenerator();
		URNGenerator(double seed);
	};

	N_Vector init_model_vector(Model &model, URNGenerator urn);

}
