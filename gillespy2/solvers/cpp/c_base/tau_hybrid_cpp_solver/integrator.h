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

	/* IntegratorStatus: represents the runtime state of the integrator.
	 * OK indicates that no errors have occurred.
	 */
	enum IntegrationStatus
	{
		// No errors have occurred.
		OK = 0,
		// Attempted to perform a SUNDIALS operation on a null CVODE object.
		NULL_POINTER,
		// A non-null object resulted in a memory error and must be initialized.
		BAD_MEMORY,
		// Could not perform integration, step size too small.
		BAD_STEP_SIZE
	};

	struct IntegratorData
	{
		HybridSimulation *simulation;
		std::vector<HybridSpecies> *species_state;
		std::vector<HybridReaction> *reaction_state;

		std::vector<double> concentrations;
		std::vector<int> populations;
		std::vector<double> propensities;

		IntegratorData(HybridSimulation *simulation);
		IntegratorData(HybridSimulation *simulation, int num_species, int num_reactions);
		IntegratorData(IntegratorData &prev_data);
	};

	/* :IntegrationResults:
	 * Organized data structure for accessing the integrator's output vector.
	 * Contents are MUTABLE! Updating the values in any containing pointers
	 *   will be permanently reflected in the integrator's vector.
	 * 
	 * All pointers in the structure point to different regions of the same vector.
	 * N_Vector: [ --- concentrations --- | ---- rxn_offsets ---- ]
	 */
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
		N_Vector y0;
		double t0;
		SUNLinearSolver solver;
		int num_species;
		int num_reactions;
	public:
		// status: check for errors before using the results.
		IntegrationStatus status;
		N_Vector y;
		realtype t;

		/* save_state()
		 * Creates a duplicate copy of the integrator's current solution vector.
		 * Contents of the most recent duplicate will be restored when restore_state() is called.
		 * 
		 * Returns the time value of the integrator's saved state.
		 */
		double save_state();

		/* restore_state()
		 * Loads the most recent duplicated copy of the solution vector.
		 * 
		 * Returns the time value that the integrator was restored to.
		 */
		double restore_state();

		/* refresh_state()
		 * Loads any new changes to the solution vector without changing previous output.
		 * Any new values assigned to the public N_Vector y will be recognized by the integrator.
		 * The current time value remains the same. To change this, modify `t`.
		 */
		void refresh_state();

		void reinitialize(N_Vector y_reset);

		IntegrationResults integrate(double *t);
		IntegratorData data;

		Integrator(HybridSimulation *simulation, N_Vector y0, double reltol, double abstol);
		~Integrator();
	};

	struct URNGenerator
	{
	private:
		std::uniform_real_distribution<double> uniform;
		std::mt19937_64 rng;
		unsigned long long seed;
	public:
		double next();
		URNGenerator();
		explicit URNGenerator(unsigned long long seed);
	};

	N_Vector init_model_vector(Model<double> &model, URNGenerator urn);

	int rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data);

}
