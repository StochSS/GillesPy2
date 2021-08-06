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

#include "model.h"
#include "tau.h"

#include <vector>

#define GPY_HYBRID_ABSTOL 1e-8
#define GPY_HYBRID_RELTOL 1e-8

namespace Gillespy
{
	namespace TauHybrid
	{
		typedef int ReactionId;

		struct EventOutput
		{
			double *species_out;
			double *variable_out;
			const double *species;
			const double *variables;
			const double *constants;
		};

		class EventExecution;

		class Event
		{
		public:
			friend class EventExecution;

			inline bool trigger(double t, const double *state) const
			{
				return Event::trigger(m_event_id, t, state, Reaction::s_variables.get(), Reaction::s_constants.get());
			}

			inline double delay(double t, const double *state) const
			{
				return Event::delay(m_event_id, t, state, Reaction::s_variables.get(), Reaction::s_constants.get());
			}

			EventExecution get_execution(double t,
					const double *state, int num_state) const;
			static void use_events(std::vector<Event> &events);

		private:
			int m_event_id;
			bool m_use_trigger_state;
			bool m_use_persist;

			explicit Event(int event_id, bool use_trigger_state, bool use_persist);

			static bool trigger(
					int event_id, double t,
					const double *state,
					const double *variables,
					const double *constants);
			static double delay(
					int event_id, double t,
					const double *state,
					const double *variables,
					const double *constants);
			static double priority(
					int event_id, double t,
					const double *state,
					const double *variables,
					const double *constants);
			static void assign(int event_id, double t, EventOutput output);
		};

		class EventExecution
		{
		public:

			friend class Event;
			~EventExecution();
			EventExecution(const EventExecution&);
			EventExecution(EventExecution&&) noexcept;
			EventExecution &operator=(const EventExecution&);
			EventExecution &operator=(EventExecution&&) noexcept;

			void execute(double t, EventOutput output) const;
			void execute(double t, double *state);
			inline double priority(double t, const double *state) const
			{
				return Event::priority(m_event_id, t, state, Reaction::s_variables.get(), Reaction::s_constants.get());
			}

			inline double get_execution_time() const { return m_execution_time; }

			bool operator<(const EventExecution &rhs) const;
			bool operator>(const EventExecution &rhs) const;

		private:
			double m_execution_time;
			int m_event_id;

			int m_num_state = 0;
			double *m_state = nullptr;

			int m_num_variables = 0;
			double *m_variables = nullptr;

			std::vector<int> m_assignments;
			void use_assignments();

			EventExecution(int event_id, double t);
			EventExecution(int event_id, double t,
						   const double *state, int num_state,
						   const double *variables, int num_variables);
		};

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
			std::vector<std::function<double(double *, int *)>> formulas;
			std::vector<std::function<double(double, double *, double *, const double *)>> rate_rules;

			double evaluate(double t, double *ode_state, int *ssa_state);
		};

		enum SimulationState : unsigned int
		{
			CONTINUOUS = 0,
			DISCRETE = 1,
			DYNAMIC = 2
		};

		struct HybridSpecies
		{
			Species<double> *base_species;

			// allows the user to specify if a species' population should definitely be modeled continuously or
			// discretely
			// CONTINUOUS or DISCRETE
			// otherwise, mode will be determined by the program (DYNAMIC)
			// if no choice is made, DYNAMIC will be assumed
			SimulationState user_mode;

			// during simulation execution, a species will fall into either of the two categories, CONTINUOUS or DISCRETE
			// this is pre-determined only if the user_mode specifies CONTINUOUS or DISCRETE.
			// otherwise, if DYNAMIC is specified, partition_mode will be continually calculated throughout the simulation
			// according to standard deviation and coefficient of variance.
			SimulationState partition_mode;

			// Tolerance level for considering a dynamic species deterministically, value is compared
			// to an estimated sd/mean population of a species after a given time step.
			//  This value will be used if a switch_min is not provided. The default value is 0.03
			double switch_tol = 0.03;

			//Minimum population value at which species will be represented as continuous.
			// If a value is given, switch_min will be used instead of switch_tol.
			unsigned int switch_min = 0;

			DifferentialEquation diff_equation;

			// Boundary condition species are not directly updated by reactions, while standard ones are.
			// If `boundary_condition` is true, then reactants are not consumed, and products are not produced.
			bool boundary_condition = false;

			HybridSpecies();
		};

		struct HybridReaction
		{
			Reaction *base_reaction;
			SimulationState mode;

			HybridReaction();
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
}
