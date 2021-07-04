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

#include <functional>
#include "model.h"
#include "tau.h"

#define GPY_HYBRID_ABSTOL 1e-5
#define GPY_HYBRID_RELTOL 1e-5

namespace Gillespy::TauHybrid
{

	typedef int ReactionId;

	/* Gillespy::TauHybrid::DiffEquation
	 * A vector containing evaluable functions, which accept integrator state and return propensities.
	 * 
	 * The vector is understood to be an arbitrarily sized collection of propensity evaluations,
	 *   each weighted by some individual, constant factor.
	 * The sum of evaulations of all collected functions is interpreted to be the dydt of that state.
	 */
	struct DifferentialEquation
	{
	public:
		std::vector<std::function<double(std::vector<double>&, std::vector<int>&)>> formulas;
		double evaluate(std::vector<double> &ode_state, std::vector<int> &ssa_state);
	};

	enum SimulationState
	{
		CONTINUOUS = 0,
		DISCRETE,
		DYNAMIC
	};

	struct HybridSpecies
	{
		Species<double> *base_species;

		// allows the user to specify if a species' population should definitely be modeled continuously or 
		// discretely
		// CONTINUOUS or DISCRETE
		// otherwise, mode will be determined by the program (DYNAMIC)
		// if no choice is made, DYNAMIC will be assumed 
		SimulationState  user_mode : 2;

		// during simulation execution, a species will fall into either of the two categories, CONTINUOUS or DISCRETE
		// this is pre-determined only if the user_mode specifies CONTINUOUS or DISCRETE.
		// otherwise, if DYNAMIC is specified, partition_mode will be continually calculated throughout the simulation
		// according to standard deviation and coefficient of variance.
		SimulationState partition_mode : 1;

		// Tolerance level for considering a dynamic species deterministically, value is compared
		// to an estimated sd/mean population of a species after a given time step.
		//  This value will be used if a switch_min is not provided. The default value is 0.03
		double switch_tol = 0.03;

		//Minimum population value at which species will be represented as continuous. 
		// If a value is given, switch_min will be used instead of switch_tol.
		unsigned int switch_min = 0;

		DifferentialEquation diff_equation;

		HybridSpecies();
	};

	struct HybridReaction
	{
		Reaction *base_reaction;
		SimulationState mode : 1;

		HybridReaction();

		double ode_propensity(ReactionId reaction_number, std::vector<double> &state);
		double ssa_propensity(ReactionId reaction_number, std::vector<int> &state);
	};

	struct HybridSimulation : Simulation<double>
	{
		std::vector<HybridSpecies> species_state;
		std::vector<HybridReaction> reaction_state;

		HybridSimulation();
		HybridSimulation(const Model<double> &model);
	};

	std::set<int> flag_det_rxns(
		std::vector<HybridReaction> &reactions,
		std::vector<HybridSpecies> &species);

	void partition_species(
		std::vector<HybridReaction> &reactions,
		std::vector<HybridSpecies> &species,
		const std::vector<double> &propensity_values,
		std::vector<double> &curr_state,
		double tau_step,
		const TauArgs<double> &TauArgs);

	void update_species_state(
		std::vector<HybridSpecies> &species,
		std::vector<double> &current_state);

	void create_differential_equations(
		std::vector<HybridSpecies> &species,
		std::vector<HybridReaction> &reactions);

}
