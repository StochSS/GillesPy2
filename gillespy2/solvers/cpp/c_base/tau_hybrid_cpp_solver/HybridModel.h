/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2022 GillesPy2 developers.
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
#include <queue>
#include <list>

namespace Gillespy
{
    namespace TauHybrid
    {
        class EventExecution;
        class Event;

        typedef int ReactionId;
        template <typename T>
        using DelayedExecutionQueue = std::priority_queue<EventExecution, std::vector<EventExecution>, T>;

        struct EventOutput
        {
            double *species_out;
            double *variable_out;
            const double *species;
            const double *variables;
            const double *constants;
        };

        class EventList
        {
        public:
            EventList();
            bool evaluate_triggers(double *event_state, double t);
            bool evaluate(double *output, int output_size, double t, const std::set<int> &events_found);
            inline bool has_active_events() const
            {
                return !m_trigger_pool.empty();
            }

        private:
            std::vector<Event> m_events;
            std::set<int> m_trigger_pool;
            std::map<int, bool> m_trigger_state;
            DelayedExecutionQueue<std::greater<EventExecution>> m_delay_queue;
            std::list<EventExecution> m_volatile_queue;
        };

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

            inline bool get_initial_value() const
            {
                return Event::initial_value(m_event_id);
            }

            inline bool is_persistent() const { return m_use_persist; }
            inline int get_event_id() const { return m_event_id; }

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
            static bool initial_value(int event_id);
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
            inline bool trigger(double t, const double *state) const
            {
                return Event::trigger(m_event_id, t, state, Reaction::s_variables.get(), Reaction::s_constants.get());
            }

            inline double get_execution_time() const { return m_execution_time; }
            inline int get_event_id() const { return m_event_id; }

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
         * The sum of evaluations of all collected functions is interpreted to be the dy/dt of that state.
         */
        struct DifferentialEquation
        {
        public:
            std::vector<std::function<double(double *)>> formulas;
            std::vector<std::function<double(double, double *, double *, const double *)>> rate_rules;

            double evaluate(double t, double *ode_state);
        };

        enum SimulationState : unsigned int
        {
            CONTINUOUS = 0,
            DISCRETE = 1,
            DYNAMIC = 2
        };

        struct HybridSpecies
        {
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

            HybridSpecies() = delete;
            explicit HybridSpecies(Species<double> *species);

        private:
            Species<double> *m_base_species;
        };

        struct HybridReaction
        {
            SimulationState mode;

            inline double propensity(double *state) const
            {
                return propensity(state, Reaction::s_variables.get(), Reaction::s_constants.get());
            }

            inline double propensity(double *state, double *parameters, const double *constants) const
            {
                switch (mode)
                {
                case SimulationState::CONTINUOUS:
                    return ode_propensity(state, parameters, constants);
                case SimulationState::DISCRETE:
                default:
                    return ssa_propensity(state, parameters, constants);
                }
            }

            inline double ode_propensity(double *state) const
            {
                return ode_propensity(state, Reaction::s_variables.get(), Reaction::s_constants.get());
            }

            inline double ode_propensity(double *state, double *parameters, const double *constants) const
            {
                return map_ode_propensity(m_id, state, parameters, constants);
            }

            inline double ssa_propensity(double *state) const
            {
                return ssa_propensity(state, Reaction::s_variables.get(), Reaction::s_constants.get());
            }

            inline double ssa_propensity(double *state, double *parameters, const double *constants) const
            {
                return map_ssa_propensity(m_id, state, parameters, constants);
            }

            inline void set_base_reaction(Reaction *reaction)
            {
                m_base_reaction = reaction;
                m_id = reaction->id;
            }

            inline Reaction *get_base_reaction() const
            {
                return m_base_reaction;
            }

            HybridReaction() = delete;
            explicit HybridReaction(Reaction *base_reaction);

        private:
            Reaction *m_base_reaction;
            ReactionId m_id;
        };

        struct HybridSimulation : Simulation<double>
        {
            std::vector<HybridSpecies> species_state;
            std::vector<HybridReaction> reaction_state;

            HybridSimulation();

            HybridSimulation(const Model<double> &model);

            enum SIMULATION_STATUS : uint8_t
            {
                UNKNOWN = 1,
                LOOP_OVER_INTEGRATE = 2,
                INTEGRATOR_FAILED = 3,
                INVALID_AFTER_SSA = 4,
                NEGATIVE_STATE_NO_SSA_REACTION = 5,
                NEGATIVE_STATE_AT_BEGINING_OF_STEP = 6,
            };
        };

        std::set<int> flag_det_rxns(
                std::vector<HybridReaction> &reactions,
                std::vector<HybridSpecies> &species);

        void partition_species(
                double current_time,
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
