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

#include "HybridModel.h"

#include <cstring>

namespace Gillespy
{
    namespace TauHybrid
    {
        Event::Event(int event_id, bool use_trigger_state, bool use_persist)
                : m_event_id(event_id),
                  m_use_trigger_state(use_trigger_state),
                  m_use_persist(use_persist)
        {}

        EventExecution Event::get_execution(double t,
                                            const double *state, int num_state) const
        {
            return m_use_trigger_state
                   ? EventExecution(m_event_id, t, state, num_state,
                                    Reaction::s_variables.get(), Reaction::s_num_variables)
                   : EventExecution(m_event_id, t);
        }


        EventExecution::EventExecution(int event_id, double t)
                : m_execution_time(t),
                  m_event_id(event_id)
        {
            use_assignments();
        }

        EventExecution::EventExecution(int event_id, double t,
           const double *state, int num_state,
           const double *variables, int num_variables)
               : m_execution_time(t),
                 m_event_id(event_id),
                 m_num_state(num_state),
                 m_state(new double[num_state]),
                 m_num_variables(num_variables),
                 m_variables(new double[num_variables])
        {
            std::memcpy(m_state, state, sizeof(double) * num_state);
            std::memcpy(m_variables, variables, sizeof(double) * num_variables);
            use_assignments();
        }

        EventExecution::EventExecution(const EventExecution &old_event)
                : m_num_state(old_event.m_num_state),
                  m_num_variables(old_event.m_num_variables),
                  m_execution_time(old_event.m_execution_time),
                  m_event_id(old_event.m_event_id),
                  m_assignments(old_event.m_assignments)
        {
            if (old_event.m_state != nullptr && m_num_state > 0)
            {
                m_state = new double[m_num_state];
                std::memcpy(m_state, old_event.m_state, sizeof(double) * m_num_state);
            }

            if (old_event.m_variables != nullptr && m_num_variables > 0)
            {
                m_variables = new double[m_num_variables];
                std::memcpy(m_variables, old_event.m_variables, sizeof(double) * m_num_variables);
            }
        }

        EventExecution::EventExecution(EventExecution &&old_event) noexcept
                : m_state(old_event.m_state),
                  m_num_state(old_event.m_num_state),
                  m_variables(old_event.m_variables),
                  m_num_variables(old_event.m_num_variables),
                  m_execution_time(old_event.m_execution_time),
                  m_event_id(old_event.m_event_id),
                  m_assignments(std::move(old_event.m_assignments))
        {
            old_event.m_num_state = 0;
            old_event.m_state = nullptr;
            old_event.m_num_variables = 0;
            old_event.m_variables = nullptr;
        }

        EventExecution &EventExecution::operator=(const EventExecution &old_event)
        {
            m_execution_time = old_event.m_execution_time;
            m_assignments = old_event.m_assignments;

            if (this != &old_event)
            {
                if (old_event.m_state != nullptr)
                {
                    // If the containers are not of equal size, then we cannot reuse heap data.
                    if (m_num_state != old_event.m_num_state)
                    {
                        delete[] m_state;
                        m_num_state = old_event.m_num_state;
                        m_state = new double[m_num_state];
                    }
                    std::memcpy(m_state, old_event.m_state, sizeof(double) * m_num_state);
                }

                if (old_event.m_variables != nullptr)
                {
                    if (m_num_variables != old_event.m_num_variables)
                    {
                        delete[] m_variables;
                        m_num_variables = old_event.m_num_variables;
                        m_variables = new double[m_num_variables];
                    }
                    std::memcpy(m_variables, old_event.m_variables, sizeof(double) * m_num_variables);
                }
            }

            return *this;
        }

        EventExecution &EventExecution::operator=(EventExecution &&old_event) noexcept
        {
            m_execution_time = old_event.m_execution_time;
            m_assignments = std::move(old_event.m_assignments);

            if (this != &old_event)
            {
                m_num_state = old_event.m_num_state;
                m_state = old_event.m_state;
                old_event.m_num_state = 0;
                old_event.m_state = nullptr;

                m_num_variables = old_event.m_num_variables;
                m_variables = old_event.m_variables;
                old_event.m_num_variables = 0;
                old_event.m_variables = nullptr;
            }

            return *this;
        }

        EventExecution::~EventExecution()
        {
            delete[] m_state;
            delete[] m_variables;
        }

        void EventExecution::execute(double t, EventOutput output) const
        {
            for (int assign_id : m_assignments)
            {
                Event::assign(assign_id, t, output);
            }
        }

        void EventExecution::execute(double t, double *state)
        {
            if (m_state == nullptr || m_variables == nullptr)
            {
                execute(t, EventOutput {
                    state,
                    Reaction::s_variables.get(),
                    state,
                    Reaction::s_variables.get(),
                    Reaction::s_constants.get()
                });
            }
            else
            {
                execute(t, EventOutput {
                    state,
                    Reaction::s_variables.get(),
                    m_state,
                    m_variables,
                    Reaction::s_constants.get()
                });
            }
        }

        bool EventExecution::operator<(const EventExecution &rhs) const
        {
            return m_execution_time < rhs.m_execution_time;
        }

        bool EventExecution::operator>(const EventExecution &rhs) const
        {
            return m_execution_time > rhs.m_execution_time;
        }

        HybridReaction::HybridReaction(Reaction *base_reaction)
            : mode(SimulationState::DISCRETE),
              m_base_reaction(base_reaction)
        {
            m_id = base_reaction->id;
        }

        HybridSpecies::HybridSpecies(Species<double> *base_species)
                : user_mode(SimulationState::DYNAMIC),
                  partition_mode(SimulationState::DISCRETE),
                  m_base_species(base_species)
        {
            // Empty constructor body
        }

        HybridSimulation::HybridSimulation()
                : Simulation<double>()
        {
            // Empty constructor body
        }

        HybridSimulation::HybridSimulation(const Model<double> &model)
                : Simulation<double>()
        {
            for (int spec_i = 0; spec_i < model.number_species; ++spec_i)
            {
                species_state.emplace_back(HybridSpecies(model.species.get() + spec_i));
            }

            for (int rxn_i = 0; rxn_i < model.number_reactions; ++rxn_i)
            {
                reaction_state.emplace_back(HybridReaction(model.reactions.get() + rxn_i));
            }
        }


        double DifferentialEquation::evaluate(
                const double t,
                double *ode_state)
        {
            double sum = 0.0;

            for (auto &rate_rule : rate_rules)
            {
                sum += rate_rule(t, ode_state, Reaction::s_variables.get(), Reaction::s_constants.get());
            }

            for (auto &formula : formulas)
            {
                sum += formula(ode_state);
            }

            return sum;
        }


        void create_differential_equations(
                std::vector<HybridSpecies> &species,
                std::vector<HybridReaction> &reactions)
        {
            // For now, differential equations are generated from scratch.
            // It may be more efficient to determine which formulas need to change.
            // Until then, the compound formulas in every species are cleared.
            for (HybridSpecies &spec : species)
            {
                spec.diff_equation.formulas.clear();
            }

            for (HybridReaction &rxn : reactions)
            {
                if (rxn.mode == SimulationState::DISCRETE)
                {
                    continue;
                }

                for (int spec_i = 0; spec_i < species.size(); ++spec_i)
                {
                    // A species change of 0 indicates that this species is not a dependency for this reaction.
                    if (rxn.get_base_reaction()->species_change[spec_i] == 0)
                    {
                        continue;
                    }

                    HybridSpecies &spec = species[spec_i];
                    auto &formula_set = spec.diff_equation.formulas;
                    int spec_diff = rxn.get_base_reaction()->species_change[spec_i];

                    formula_set.emplace_back([&rxn, spec_diff](double *state) {
                        return spec_diff * rxn.propensity(state);
                    });
                }
            }
        }

        // Helper method to flag reactions that can be processed deterministically (continuous change)
        // without exceeding the user-supplied tolerance
        std::set<int> flag_det_rxns(
                std::vector<HybridReaction> &reactions,
                std::vector<HybridSpecies> &species)
        {
            int num_reactions = reactions.size();
            int num_species = species.size();
            std::set<int> det_rxns;

            for (int rxn_i = 0; rxn_i < num_reactions; ++rxn_i)
            {
                // start with the assumption that reaction is determinstic
                HybridReaction &rxn = reactions[rxn_i];
                rxn.mode = SimulationState::CONTINUOUS;

                // iterate through the dependent species of this reaction
                for (int spec_i = 0; spec_i < num_species && rxn.mode == SimulationState::CONTINUOUS; ++spec_i)
                {
                    if (rxn.get_base_reaction()->reactants_change[spec_i] == 0 && 
                        rxn.get_base_reaction()->products_change[spec_i] == 0) {
                        continue;
                    }

                    // if any of the dependencies are set by the user as discrete OR
                    // have been set as dynamic and has not been flagged as deterministic,
                    // allow it to be modelled discretely
                    if (species[spec_i].user_mode == SimulationState::DYNAMIC) {
                        rxn.mode = species[spec_i].partition_mode;
                    } else {
                        rxn.mode = species[spec_i].user_mode;
                    }
                    // Loop breaks if we've already determined that it is to be marked as discrete.
                    if(rxn.mode == SimulationState::DISCRETE){
                        break;
                    }
                }

                if (rxn.mode == SimulationState::CONTINUOUS) {
                    det_rxns.insert(rxn_i);
                }
            }

            return det_rxns;
        }

        void partition_species(
                double current_time,
                std::vector<HybridReaction> &reactions,
                std::vector<HybridSpecies> &species,
                const std::vector<double> &propensity_values,
                std::vector<double> &curr_state,
                double tau_step,
                const TauArgs<double> &tauArgs)
        {
            

            // coefficient of variance- key:species id, value: cv
            std::map<int, double> CV;
            // means
            std::map<int, double> means;
            // standard deviation
            std::map<int, double> sd;

            // Initialize means and sd's
            for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
                HybridSpecies &spec = species[spec_i];

                if (spec.user_mode == SimulationState::DYNAMIC) {
                    means.insert({spec_i, curr_state[spec_i]});
                    sd.insert({spec_i, 0});
                }
            }

            // calculate means and standard deviations for dynamic-mode species involved in reactions
            for (int rxn_i = 0; rxn_i < reactions.size(); ++rxn_i) {
                HybridReaction &rxn = reactions[rxn_i];
                for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
                    HybridSpecies &spec = species[spec_i];
                    // Only dynamic species are to be considered.
                    if (spec.user_mode != SimulationState::DYNAMIC) {
                        continue;
                    }
                    // Selected species is either a reactant or a product, 
                    if(rxn.get_base_reaction()->reactants_change[spec_i] > 0) {
                        // Selected species is a reactant.
                        means[spec_i] -= (propensity_values[rxn_i] * rxn.get_base_reaction()->reactants_change[spec_i]);
                        sd[spec_i] += (propensity_values[rxn_i] * std::pow(rxn.get_base_reaction()->reactants_change[spec_i], 2));
                    } 
                    if(rxn.get_base_reaction()->products_change[spec_i] > 0) {
                        // Selected species is a product.
                        HybridSpecies &product = species[spec_i];
                        means[spec_i] += (propensity_values[rxn_i] * rxn.get_base_reaction()->products_change[spec_i]);
                        sd[spec_i] += (propensity_values[rxn_i] * std::pow(rxn.get_base_reaction()->products_change[spec_i], 2));
                    }
                }
            }
            // calculate coefficient of variation using means and sd
            for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
                HybridSpecies &spec = species[spec_i];
                // Only dynamic species are to be considered.
                if (spec.user_mode != SimulationState::DYNAMIC) {
                    continue;
                }
                if (means[spec_i] > 0 && sd[spec_i] > 0) {
                    CV[spec_i] = (sqrt(sd[spec_i]) / means[spec_i]);
                } else {
                    CV[spec_i] = 1;
                }
            }

            // keep a history of the past CV values, and calculate a time-average values
            std::map<int, double> CV_a;
            static std::map<int, std::queue<double>> cv_history;
            static std::map<int, double> cv_history_sum;
            int history_length = 12;
            // 
            if( current_time == 0.0 ){  // reset cv_history at start of trajectory
                cv_history.clear();
                cv_history_sum.clear();
            }

            for (int spec_i = 0; spec_i < species.size(); ++spec_i)
            {
                
                HybridSpecies &spec = species[spec_i];
                // Only dynamic species are to be considered.
                if (spec.user_mode != SimulationState::DYNAMIC){
                    continue;
                }
                if(cv_history.count(spec_i) == 0){
                    cv_history[spec_i] = std::queue<double>();
                }
                if(cv_history_sum.count(spec_i) == 0){
                    cv_history_sum[spec_i] = 0.0;
                }
                cv_history[spec_i].push(CV[spec_i]);
                cv_history_sum[spec_i] += CV[spec_i];
                if(cv_history[spec_i].size() > history_length){
                    double removed = cv_history[spec_i].front();
                    cv_history[spec_i].pop();
                    cv_history_sum[spec_i] -= removed;
                }
                CV_a[spec_i] = cv_history_sum[spec_i]/cv_history[spec_i].size();
            }

            // Select DISCRETE or CONTINOUS mode for each species
            for (int spec_i = 0; spec_i < species.size(); ++spec_i) {
                HybridSpecies &spec = species[spec_i];

                // Only dynamic species are to be considered.
                if (spec.user_mode != SimulationState::DYNAMIC) {
                    continue;
                }
                unsigned int prev_partition = spec.partition_mode;
                if (spec.switch_min == 0) {
                    // (default value means switch  min not set, use switch tol)
                    spec.partition_mode = CV_a[spec_i] < spec.switch_tol
                                          ? SimulationState::CONTINUOUS
                                          : SimulationState::DISCRETE;
                } else {
                    spec.partition_mode = means[spec_i] > spec.switch_min
                                          ? SimulationState::CONTINUOUS
                                          : SimulationState::DISCRETE;
                }
            }
        }

        void update_species_state(
                std::vector<HybridSpecies> &species,
                std::vector<double> &current_state)
        {
            for (int spec_i = 0; spec_i < species.size(); ++spec_i)
            {
                switch (species[spec_i].partition_mode)
                {
                case SimulationState::DISCRETE:
                    current_state[spec_i] = std::round(current_state[spec_i]);
                    break;
                default:
                    break;
                }
            }
        }

        EventList::EventList()
        {
            Event::use_events(m_events);

            for (auto &event : m_events)
            {
                // With the below implementation, it is impossible for an event to fire at t=t[0].
                m_trigger_state.insert({
                    event.get_event_id(),
                    event.get_initial_value(),
                });
            }
        }

        bool EventList::evaluate_triggers(double *event_state, double t)
        {
            for (auto &event: m_events)
            {
                if (event.trigger(t, event_state) != m_trigger_state.at(event.get_event_id()))
                {
                    m_trigger_pool.insert(event.get_event_id());
                }
            }

            return has_active_events();
        }

        bool EventList::evaluate(double *event_state, int output_size, double t, const std::set<int> &events_found)
        {
            if (m_events.empty())
            {
                return has_active_events();
            }

            auto compare = [t, event_state](EventExecution &lhs, EventExecution &rhs) -> bool
            {
                return lhs.priority(t, event_state) < rhs.priority(t, event_state);
            };
            std::priority_queue<EventExecution, std::vector<EventExecution>, decltype(compare)>
                    trigger_queue(compare);

            // Step 1: Identify any fired event triggers.
            for (auto &event : m_events)
            {
                if (m_trigger_state.at(event.get_event_id()) != event.trigger(t, event_state))
                {
                    double delay = event.delay(t, event_state);
                    EventExecution execution = event.get_execution(t + delay, event_state, output_size);

                    if (delay <= 0)
                    {
                        // Immediately put EventExecution on "triggered" pile
                        trigger_queue.push(execution);
                    }
                    else if (event.is_persistent())
                    {
                        // Put EventExecution on "delayed" pile
                        m_delay_queue.push(execution);
                    }
                    else
                    {
                        // Search the volatile queue to see if it is already present.
                        // If it is, the event has "double-fired" and must be erased.
                        auto vol_iter = m_volatile_queue.begin();
                        while (vol_iter != m_volatile_queue.end()
                               && vol_iter->get_event_id() != event.get_event_id())
                        {
                            ++vol_iter;
                        }

                        if (vol_iter == m_volatile_queue.end())
                        {
                            // No match found; this is a new delay trigger, and is therefore valid.
                            // Delayed, but must be re-checked on every iteration.
                            m_volatile_queue.push_back(execution);
                        }
                        else
                        {
                            // Match found; this is an existing trigger, discard.
                            m_volatile_queue.erase(vol_iter);
                            m_trigger_pool.erase(event.get_event_id());
                            m_trigger_state.at(event.get_event_id()) = !m_trigger_state.at(event.get_event_id());
                        }
                    }
                }
            }

            // Step 2: Process delayed, non-persistent executions that are now ready to fire.
            // Both the volatile and non-volatile queue are processed in a similar manner.
            for (auto vol_event = m_volatile_queue.begin(); vol_event != m_volatile_queue.end(); ++vol_event)
            {
                // Execution objects in the volatile queue must remain True until execution.
                // Remove any execution objects which transitioned to False before execution.
                if (vol_event->get_execution_time() < t)
                {
                    trigger_queue.push(*vol_event);
                    vol_event = m_volatile_queue.erase(vol_event);
                }
            }

            // Step 3: Process delayed executions, which includes both persistent triggers
            // and non-persistent triggers whose execution time has arrived.
            while (!m_delay_queue.empty())
            {
                auto &event = m_delay_queue.top();
                if (event.get_execution_time() >= t)
                {
                    // Delay queue is sorted in chronological order.
                    // As soon as we hit a time that is beyond the current time,
                    //  there is no use in continuing through the queue.
                    break;
                }
                trigger_queue.push(event);
                m_delay_queue.pop();
            }

            // Step 4: Process any pending triggers, unconditionally.
            while (!trigger_queue.empty())
            {
                auto event = trigger_queue.top();

                event.execute(t, event_state);
                trigger_queue.pop();
                m_trigger_pool.erase(event.get_event_id());
            }

            // Step 5: Update trigger states based on the new simulation state.
            // This is to account for events that re-assign values that the event triggers depend on.
            for (auto &event : m_events)
            {
                m_trigger_state.at(event.get_event_id()) = event.trigger(t, event_state);
            }

            return has_active_events();
        }
    }
}
