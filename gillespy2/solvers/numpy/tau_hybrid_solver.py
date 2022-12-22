# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import copy
import random, math, sys
from collections import OrderedDict
from scipy.integrate import ode, LSODA
import heapq
import numpy as np
import threading
import gillespy2
from gillespy2.solvers.utilities import Tau
from gillespy2.core import GillesPySolver, log, Event, RateRule, AssignmentRule, FunctionDefinition
from gillespy2.core.gillespyError import *
from gillespy2.core.results import Results

eval_globals = math.__dict__

def __piecewise(*args):
    # Eval entry for piecewise functions
    args = list(args)
    sol = None
    if len(args) % 2:
        args.append(True)
    for i, arg in enumerate(args):
        if not i % 2:
            continue
        if arg:
            sol = args[i - 1]
            break
    return sol


def __xor(*args):
    # Eval entry for MathML xor function
    from operator import ixor
    from functools import reduce
    args = list(args)
    return reduce(ixor, args)


eval_globals['false'] = False
eval_globals['true'] = True
eval_globals['piecewise'] = __piecewise
eval_globals['xor'] = __xor


class TauHybridSolver(GillesPySolver):
    """
    This solver uses a root-finding interpretation of the direct SSA method,
    along with ODE solvers to simulate ODE and Stochastic systems
    interchangeably or simultaneously.
    Uses integrators from scipy.integrate.ode to perform calculations used to produce solutions.

    :param model: The model on which the solver will operate.
    :type model: gillespy2.Model
    """
    name = "TauHybridSolver"
    rc = 0
    result = None
    stop_event = None

    def __init__(self, model=None, profile_reactions=False):
        if model is None:
            raise SimulationError("A model is required to run the simulation.")

        name = 'TauHybridSolver'
        rc = 0
        self.model = copy.deepcopy(model)
        self.is_instantiated = True
        self.profile_reactions = profile_reactions
        self.profile_data = {}
        if self.profile_reactions:
            self.profile_data['time'] = []
            for k in list(self.model.listOfSpecies)+list(self.model.listOfReactions):
                self.profile_data[k] = []
        self.non_negative_species = set()
        for reaction, _ in model.listOfReactions.items():
            for key, value in model.listOfReactions[reaction].reactants.items():
                self.non_negative_species.add(key.name)
            for key, value in model.listOfReactions[reaction].products.items():
                self.non_negative_species.add(key.name)

    def __save_state_to_output(self, curr_time, save_index, curr_state, species, 
                                trajectory, save_times):
        """
        Helper function to save the curr_state to the trajectory output
        """
        
        # Now update the step and trajectories for this step of the simulation.
        # Here we make our final assignments for this step, and begin
        # populating our results trajectory.
           
        num_saves = 0
        for time in save_times:
            if time > curr_time:
                break
            # if a solution is given for it
            trajectory_index = save_index
            assignment_state = copy.deepcopy(curr_state)
            for s,sname in enumerate(species):
                # Get ODE Solutions
                trajectory[trajectory_index][s + 1] = curr_state[sname]
                # Update Assignment Rules for all processed time points
                if len(self.model.listOfAssignmentRules):
                    # Copy ODE state for assignments
                    assignment_state[sname] = curr_state[sname]
            assignment_state['t'] = time
            for ar in self.model.listOfAssignmentRules.values():
                assignment_value = eval(ar.formula, {**eval_globals, **assignment_state})
                assignment_state[ar.variable] = assignment_value
                trajectory[trajectory_index][species.index(ar.variable.name) + 1] = assignment_value
            num_saves += 1
            save_index += 1
        save_times = save_times[num_saves:]  # remove completed save times
        return save_times, save_index


    def __toggle_reactions(self, all_compiled, deterministic_reactions, dependencies,
                            curr_state, det_spec, rr_sets):
        """
        Helper method which is used to convert reaction channels into
        rate rules, and rate rules into reaction channels, as they are switched
        dynamically throughout the simulation based on user-supplied tolerance.
        """

        # initialize variables
        inactive_reactions = all_compiled['inactive_rxns']
        rate_rules = all_compiled['rules']
        rxns = all_compiled['rxns']

        # If the set has changed, reactivate non-determinsitic reactions
        reactivate = []
        for r in inactive_reactions:
            if not r in deterministic_reactions:
                reactivate.append(r)
        for r in reactivate:
            rxns[r] = inactive_reactions.pop(r, None)

        # Deactivate Determinsitic Reactions
        for r in deterministic_reactions:
            if not r in inactive_reactions:
                inactive_reactions[r] = rxns.pop(r, None)

        # floor non-det species
        for s,d in det_spec.items():
            if not d and isinstance(curr_state[s], float):
                curr_state[s] = round(curr_state[s])

        # Check if this reaction set is already compiled and in use:
        if deterministic_reactions in rr_sets.keys():
            return rr_sets[deterministic_reactions]
        else:
        # Otherwise, this is a new determinstic reaction set that must be compiled
            tmp =  self.__create_diff_eqs(deterministic_reactions,
                                            dependencies, rr_sets, rate_rules)
            return tmp

    def __create_diff_eqs(self, comb, dependencies, rr_sets, rate_rules):
        """
        Helper method used to convert stochastic reaction descriptions into
        differential equations, used dynamically throught the simulation.
        """
        diff_eqs = OrderedDict()
        output_rules = copy.deepcopy(rate_rules)

        # Initialize sample dict
        rr_vars = {}
        for n, rr in self.model.listOfRateRules.items():
            rr_vars[rr.variable.name] = n
        for spec in self.model.listOfSpecies:
            if spec in rr_vars.keys():
                diff_eqs[self.model.listOfSpecies[spec]] = self.model.listOfRateRules[rr_vars[spec]].formula
            else:
                diff_eqs[self.model.listOfSpecies[spec]] = '0'

        # loop through each det reaction and concatenate it's diff eq for each species
        for reaction in comb:
            factor = {dep: 0 for dep in dependencies[reaction]}

            for key, value in self.model.listOfReactions[reaction].reactants.items():
                if not key.constant and not key.boundary_condition:
                    factor[key.name] -= value
            for key, value in self.model.listOfReactions[reaction].products.items():
                if not key.constant and not key.boundary_condition:
                    factor[key.name] += value

            for dep in dependencies[reaction]:
                if factor[dep] != 0:
                    if self.model.listOfSpecies[dep].mode == 'continuous':
                        diff_eqs[self.model.listOfSpecies[dep]] += ' + {0}*({1})'.format(factor[dep],
                                                           self.model.listOfReactions[reaction].ode_propensity_function)
                    else:
                        diff_eqs[self.model.listOfSpecies[dep]] += ' + {0}*({1})'.format(factor[dep],
                                                           self.model.listOfReactions[reaction].propensity_function)

        for spec in self.model.listOfSpecies:
            if diff_eqs[self.model.listOfSpecies[spec]] == '0':
                del diff_eqs[self.model.listOfSpecies[spec]]
        # create a dictionary of compiled gillespy2 rate rules
        for spec, rate in diff_eqs.items():
            output_rules[spec] = compile(gillespy2.RateRule(spec, rate).formula, '<string>', 'eval')
        rr_sets[comb] = rate_rules # save values
        return output_rules

    def __flag_det_reactions(self, det_spec, det_rxn, dependencies):
        """
        Helper method used to flag reactions that can be processed
        deterministically without exceeding the user-supplied tolerance.
        """
        # Determine if each rxn would be deterministic apart from other reactions
        prev_state = det_rxn.copy()
        for rxn in self.model.listOfReactions:
            # assume it is deterministic
            det_rxn[rxn] = True
            # iterate through the dependent species of this reaction
            for species in dependencies[rxn]:
                # if any of the dependencies are discrete or (dynamic AND the
                # species itself has not been flagged as deterministic)
                # then allow it to be modelled discretely
                if self.model.listOfSpecies[species].mode == 'discrete':
                    det_rxn[rxn] = False
                    break
                if self.model.listOfSpecies[species].mode == 'dynamic' and det_spec[species] == False:
                    det_rxn[rxn] = False
                    break

        # Create a hashable frozenset of all determinstic reactions
        deterministic_reactions = set()
        for rxn in det_rxn:
            if det_rxn[rxn]:
                deterministic_reactions.add(rxn)
        deterministic_reactions = frozenset(deterministic_reactions)
        return deterministic_reactions

    def __calculate_statistics(self, curr_time, propensities, curr_state, tau_step, det_spec, cv_history={}):
        """
        Calculates Mean, Standard Deviation, and Coefficient of Variance for each
        dynamic species, then set if species can be represented determistically

        NOTE: the argument cv_history should not be passed in, this is modified by the function
        to keep a persistent data set.
        """
        #move configuration to solver init (issue #905)
        history_length = 12 #cite: "Statistical rules of thumb" G Van Belle

        if curr_time==0.0: #re-set cv_history
            for k in list(cv_history.keys()):
                del cv_history[k]

        CV = OrderedDict()
        # calculate CV by estimating the next step
        mn = {species: curr_state[species] for species, value in
              self.model.listOfSpecies.items() if value.mode == 'dynamic'}
        sd = {species: 0 for species, value in
              self.model.listOfSpecies.items() if value.mode == 'dynamic'}

        for r, rxn in self.model.listOfReactions.items():
            for reactant in rxn.reactants:
                if reactant.mode == 'dynamic':
                    mn[reactant.name] -= (propensities[r] * rxn.reactants[reactant])
                    sd[reactant.name] += (propensities[r] * (rxn.reactants[reactant] ** 2))
            for product in rxn.products:
                if product.mode == 'dynamic':
                    mn[product.name] += (propensities[r] * rxn.products[product])
                    sd[product.name] += (propensities[r] * (rxn.products[product] ** 2))
        # Calcuate the derivative based CV
        for species,value in self.model.listOfSpecies.items():
            if value.mode == 'dynamic':
                if mn[species] > 0 and sd[species] > 0:
                    CV[species] = math.sqrt(sd[species]) / mn[species]
                else:
                    CV[species] = 1  # value chosen to guarantee species will be discrete

        # Keep a history of the past CV values, calculate a time-averaged value
        CV_a = OrderedDict() # time-averaged CV (forward derivative based)
        for species,value in self.model.listOfSpecies.items():
            if value.mode == 'dynamic':
                if species not in cv_history:
                    cv_history[species] = []
                cv_history[species].append(CV[species])
                if len(cv_history[species]) > history_length:
                    cv_history[species].pop(0) #remove the first item
                CV_a[species] = sum(cv_history[species])/len(cv_history[species])

        # Select DISCRETE or CONTINOUS mode for each species
        for species in mn:
            prev_det = det_spec[species]
            sref = self.model.listOfSpecies[species]
            if sref.switch_min == 0:
                # Set species to deterministic if CV is less than threshhold
                det_spec[species] = CV_a[species] < sref.switch_tol
            else:
                det_spec[species] = mn[species] > sref.switch_min

        return mn, sd, CV

    @staticmethod
    def __f(t, y, curr_state, species, reactions, rate_rules, propensities,
            y_map, compiled_reactions, active_rr, events, assignment_rules):
        """
        Evaluate the propensities for the reactions and the RHS of the Reactions and RateRules.
        Also evaluates boolean value of event triggers.
        """
        state_change = [0] * len(y_map)
        curr_state['t'] = t
        curr_state['time'] = t
        for item, index in y_map.items():
            if item in assignment_rules:
                curr_state[assignment_rules[item].variable] = eval(assignment_rules[item].formula,
                                                                   {**eval_globals, **curr_state})
            else:
                curr_state[item] = y[index]
        for s, rr in active_rr.items():
            try:
                state_change[y_map[s.name]] += eval(rr, {**eval_globals, **curr_state})
            except ValueError as e:
                pass
        for i, r in enumerate(compiled_reactions):
            propensities[r] = eval(compiled_reactions[r], {**eval_globals, **curr_state})
            state_change[y_map[r]] += propensities[r]
        for event in events:
            triggered = eval(event.trigger.expression, {**eval_globals, **curr_state})
            if triggered:
                state_change[y_map[event]] = 1
        return state_change

    def __find_event_time(self, sol, start, end, index, depth):
        """
        Helper method providing binary search implementation for locating
        precise event times.
        """
        dense_range = np.linspace(start, end, 3)
        mid = dense_range[1]
        if start >= mid or mid >= end or depth == 20: return end
        solutions = np.diff(sol.sol(dense_range)[-len(self.model.listOfEvents) + index])
        bool_res = [x > 0 for x in solutions]

        if bool_res[0]:  # event before mid
            depth += 1
            return self.__find_event_time(sol, dense_range[0],
                                          dense_range[1], index, depth)
        else:  # event after mid
            depth += 1
            return self.__find_event_time(sol, dense_range[1],
                                          dense_range[2], index, depth)

    def __detect_events(self, event_sensitivity, sol, delayed_events,
                        trigger_states, curr_time, curr_state):
        """
        Helper method to locate precise time of event firing.  This method
        first searches for any instance of an event using event_sensitivity to
        determine the granularity of the search.  If an event is detected, a
        binary search is then used to locate the precise time of that event.
        """
        event_times = {}
        dense_range = np.linspace(sol.t[0], sol.t[-1], len(sol.t) * event_sensitivity)
        solutions = np.diff(sol.sol(dense_range))
        for i, e in enumerate(self.model.listOfEvents.values()):
            bool_res = [x > 0 for x in solutions[i - len(self.model.listOfEvents)]]
            curr_state[e.name] = bool_res[-1]
            # Search for changes from False to True in event, record first time
            for y in range(1, len(dense_range) - 1):
                # Check Persistent Delays.  If an event is not designated as
                # persistent, and the trigger expression fails to evaluate as
                # true before assignment is carried out, remove the event from
                # the queue.
                if e.name in trigger_states and not e.trigger.persistent:
                    if not bool_res[y]:
                        delayed_events[i] = delayed_events[-1]  # move to end
                        heapq.heappop(delayed_events)
                        del trigger_states[e.name]
                        curr_state[e.name] = False
                # IF triggered from false to true, refine search
                elif bool_res[y] and dense_range[y] != curr_time and bool_res[y - 1] == 0:
                    event_time = self.__find_event_time(sol, dense_range[y - 1],
                                                        dense_range[y + 1], i, 0)
                    if event_time in event_times:
                        event_times[event_time].append(e)
                    else:
                        event_times[event_time] = [e]
                    break
        return event_times

    def __get_next_step(self, event_times, reaction_times, delayed_events,
                        sim_end, next_tau):
        """
        Helper method to determine the next action to take during simulation,
        and returns that action along with the time that it occurs.
        """

        next_event_trigger = sim_end + 2
        next_delayed_event = sim_end + 3
        curr_time = sim_end

        # event triggers
        if len(event_times):
            next_event_trigger = min(event_times)

        # delayed events
        if len(delayed_events):
            next_delayed_event = delayed_events[0][0]

        # Execute rxns/events
        next_step = {sim_end: 'end', next_tau: 'tau', next_event_trigger: 'trigger', next_delayed_event: 'delay'}

        # Set time to next action
        curr_time = min(sim_end, next_tau, next_event_trigger,
                        next_delayed_event)
        return next_step[curr_time], curr_time

    def __process_queued_events(self, event_queue, trigger_states,
                                curr_state, det_spec):
        """
        Helper method which processes the events queue. Method is primarily for
        evaluating assignments at the designated state (trigger time or event
        time).
        """
        # Process all queued events
        events_processed = []
        pre_assignment_state = curr_state.copy()
        while event_queue:
            # Get events in priority order
            fired_event = self.model.listOfEvents[heapq.heappop(event_queue)[1]]
            events_processed.append(fired_event)
            if fired_event.name in trigger_states:
                assignment_state = trigger_states[fired_event.name]
                del trigger_states[fired_event.name]
            else:
                assignment_state = pre_assignment_state
            for a in fired_event.assignments:
                # Get assignment value
                assign_value = eval(a.expression, eval_globals, assignment_state)
                # Update state of assignment variable
                curr_state[a.variable.name] = assign_value

        return events_processed

    def __handle_event(self, event, curr_state, curr_time, event_queue,
                       trigger_states, delayed_events):
        """
        Helper method providing logic for updating states based on assignments.
        """
        # Fire trigger time events immediately
        if event.delay is None:
            heapq.heappush(event_queue, (eval(event.priority), event.name))
        # Queue delayed events
        else:
            curr_state['t'] = curr_time
            curr_state['time'] = curr_time
            execution_time = curr_time + eval(event.delay, {**eval_globals, **curr_state})
            curr_state[event.name] = True
            heapq.heappush(delayed_events, (execution_time, event.name))
            if event.use_values_from_trigger_time:
                trigger_states[event.name] = curr_state.copy()
            else:
                trigger_states[event.name] = curr_state

    def __check_t0_events(self, initial_state):
        """
        Helper method for firing events who reach a trigger condition at start
        of simulation, time == 0.
        """
        # Check Event State at t==0
        species_modified_by_events = []
        t0_delayed_events = {}
        for e in self.model.listOfEvents.values():
            if not e.trigger.value:
                t0_firing = eval(e.trigger.expression, {**eval_globals, **initial_state})
                if t0_firing:
                    if e.delay is None:
                        for a in e.assignments:
                            initial_state[a.variable.name] = eval(a.expression, {**eval_globals, **initial_state})
                            species_modified_by_events.append(a.variable.name)
                    else:
                        execution_time = eval(e.delay, {**eval_globals, **initial_state})
                        t0_delayed_events[e.name] = execution_time
        return t0_delayed_events, species_modified_by_events

    def __update_stochastic_rxn_states(self, compiled_reactions, curr_state, only_update=None):
        """
        Helper method for updating the state of stochastic reactions.

        if 'only_update' is set to a reaction name, it will only reset that reaction, and only one firing
        """

        rxn_count = OrderedDict()
        species_modified = OrderedDict()
        # Update stochastic reactions

        for rxn in compiled_reactions:
            rxn_count[rxn] = 0
            if only_update is not None:  # for a single SSA step
                if rxn == only_update:
                    curr_state[rxn] = math.log(random.uniform(0, 1)) #set this value, needs to be <0
                    rxn_count[rxn] = 1
            else:  # for a normal Tau-step
                while curr_state[rxn] > 0:
                    rxn_count[rxn] += 1
                    curr_state[rxn] += math.log(random.uniform(0, 1))
            if rxn_count[rxn]:
                for reactant in self.model.listOfReactions[rxn].reactants:
                    species_modified[reactant.name] = True
                    curr_state[reactant.name] -= self.model.listOfReactions[rxn].reactants[reactant] * rxn_count[rxn]
                for product in self.model.listOfReactions[rxn].products:
                    species_modified[product.name] = True
                    curr_state[product.name] += self.model.listOfReactions[rxn].products[product] * rxn_count[rxn]
        return species_modified, rxn_count

    def __integrate(self, integrator_options, curr_state, y0, curr_time,
                    propensities, y_map, compiled_reactions,
                    active_rr, event_queue,
                    delayed_events, trigger_states,
                    event_sensitivity, tau_step ):
        """
        Helper function to perform the ODE integration of one step.  This
        method uses scipy.integrate.LSODA to get simulation data, and
        determines the next stopping point of the simulation. The state is
        updated and returned to __simulate along with curr_time and the
        solution object.
        """
        max_step_size = self.model.tspan[1] - self.model.tspan[0] / 100

        from functools import partial
        events = self.model.listOfEvents.values()
        dense_output = False
        int_args = [curr_state, self.model.listOfSpecies, self.model.listOfReactions,
                    self.model.listOfRateRules,
                    propensities, y_map,
                    compiled_reactions,
                    active_rr,
                    events,
                    self.model.listOfAssignmentRules]
        rhs = lambda t, y: TauHybridSolver.__f(t, y, *int_args)
        if 'min_step' in integrator_options:
            tau_step = max(integrator_options['min_step'], tau_step)
        else:
            tau_step = max(1e-6, tau_step)
        next_tau = curr_time + tau_step
        curr_state['t'] = curr_time
        curr_state['time'] = curr_time

        # Integrate until end or tau is reached
        loop_count = 0
        sol = LSODA(rhs, curr_time, y0, next_tau)
        counter = 0
        while sol.t < next_tau:
            counter += 1
            sol.step()


        # Update states of all species based on changes made to species through
        # ODE processes.  This will update all species whose mode is set to
        # 'continuous', as well as 'dynamic' mode species which have been
        # flagged as deterministic.
       
        for spec_name, species in self.model.listOfSpecies.items():
            if not species.constant:
                curr_state[spec_name] = sol.y[y_map[spec_name]]

        # Search for precise event times
        '''
        if len(model.listOfEvents):
            event_times = self.__detect_events(event_sensitivity, sol, delayed_events,
                                               trigger_states, curr_time, curr_state)
        else:
            event_times = {}

        # Get next tau time
        reaction_times = []
        '''
        # Set curr time to next time a change occurs in the system outside of
        # the standard ODE process.  Determine what kind of change this is,
        # and set the curr_time of simulation to restart simulation after
        # making the appropriate state changes.
        event_times = {}
        reaction_times = []
        next_step, curr_time = self.__get_next_step(event_times, reaction_times,
                                                    delayed_events,
                                                    self.model.tspan[-1], next_tau)
        curr_state['t'] = curr_time

        # Stochastic Reactions are also fired through a root-finding method
        # which mirrors the standard SSA probability.  Since we are using
        # Tau-Leaping, rather than firing reactions based on a direct poisson
        # distribution from the propensity, we use that same poisson
        # distribution to select a random negative offset for the reaction
        # state, then integrate forward along with any deterministic
        # reactions/species for a designated time (tau).  Root crossings
        # satisfy the SSA firing process, and because multiple reactions can be
        # fired in a single Tau step, a new random number will be generated and
        # added (representing a single-firing each time) until the state value
        # of the reaction is once again negative.
        for rxn in compiled_reactions:
            curr_state[rxn] = sol.y[y_map[rxn]]

        # In the case that a major change occurs to the system (outside of the
        # standard ODE process) is caused by an event trigger, we then examine
        # the event which was triggered.  if Event assignments are carried out
        # immediately, we add them to a heap in Event Priority order to be
        # immediately dealt with.  If an Event contains a delay, we can now set
        # the time of event assignment execution by evaluating the event delay
        # in the current state.
        if next_step == 'trigger':
            for event in event_times[curr_time]:
                self.__handle_event(event, curr_state, curr_time,
                                    event_queue, trigger_states, delayed_events)
        # In the case that a major change occurs to the system (outside of the
        # standard ODE process) is caused by a delayed event which has now
        # reached the designated time of execution, we add all currently
        # designated event assignments to a heap to be processed in Event
        # Priority order.
        elif next_step == 'delay':
            event = heapq.heappop(delayed_events)
            heapq.heappush(event_queue, (eval(self.model.listOfEvents[event[1]].priority), event[1]))

        return sol, curr_time


    def __simulate_negative_state_check(self, curr_state):
            # check each species to see if they are negative
            for s in self.non_negative_species:
                if curr_state[s] < 0:
                    raise SimulationError(f"Negative State detected at begining of step." \
                    " Species involved in reactions can not be negative.")

    def __simulate_invalid_state_check(self, species_modified, curr_state, compiled_reactions):
        invalid_state = False
        err_message=""
        # check each species to see if they are negative
        for s in species_modified.keys():
            if curr_state[s] < 0:
                invalid_state = True
                err_message += f"'{s}' has negative state '{curr_state[s]}'"
        return (invalid_state, err_message)

    def __simulate(self, integrator_options, curr_state, y0, curr_time,
                   propensities, species, parameters, compiled_reactions,
                   active_rr, y_map, trajectory, save_times, save_index,
                   delayed_events, trigger_states, event_sensitivity,
                   tau_step, debug, det_spec):
        """
        Function to process simulation until next step, which can be a
        stochastic reaction firing, an event trigger or assignment, or end of
        simulation.

        :param curr_state: Contains all state variables for system at current time
        :type curr_state: dict

        :param curr_time: Represents current time
        :type curr_time: float

        :param save_times: Currently unreached save points
        :type save_times: list

        :returns: curr_state, curr_time, save_times, sol
            sol - Python object returned from LSODA which contains all solution
            data.
        """

        # first check if we have a valid state:
        self.__simulate_negative_state_check(curr_state)
        if curr_time == 0.0:
            # save state at beginning of simulation
            save_times, save_index = self.__save_state_to_output(
                curr_time, save_index, curr_state, species, trajectory, save_times
            )

        event_queue = []
        prev_y0 = copy.deepcopy(y0)
        prev_curr_state = copy.deepcopy(curr_state)
        prev_curr_time = curr_time
        loop_count = 0
        invalid_state = False


        starting_curr_state = copy.deepcopy(curr_state)
        starting_propensities = copy.deepcopy(propensities)
        starting_tau_step=tau_step
        species_modified=None
        rxn_count=None
        loop_err_message=""
        if self.profile_reactions:
            self.profile_data['time'].append(curr_time)
            for k in list(self.model.listOfSpecies)+list(self.model.listOfReactions):
                self.profile_data[k].append(curr_state[k])


        # check each reaction to see if it is >=0. If we have taken a single SSA step, this could be >0 for the non-selected reactions, check if propensity is zero and reset if so
        for r in compiled_reactions.keys():
            if curr_state[r] >= 0 and propensities[r] == 0:
                curr_state[r] = math.log(random.uniform(0, 1))

        sol, curr_time = self.__integrate(integrator_options, curr_state,
                                          y0, curr_time, propensities, y_map,
                                          compiled_reactions,
                                          active_rr,
                                          event_queue,
                                          delayed_events,
                                          trigger_states,
                                          event_sensitivity,
                                          tau_step
                                          )

        species_modified,rxn_count = self.__update_stochastic_rxn_states(compiled_reactions, curr_state)

        # Occasionally, a tau step can result in an overly-aggressive
        # forward step and cause a species population to fall below 0,
        # which would result in an erroneous simulation.
        # We estimate the time to the first
        # stochatic reaction firing (assume constant propensities) and
        # simulate the ODE system until that time, fire that reaction
        # and continue the simulation.
        (invalid_state, invalid_err_message) = self.__simulate_invalid_state_check(species_modified, curr_state, compiled_reactions)

        if invalid_state:
            invalid_state = False
            # Redo this step, with a smaller Tau such that a single SSA reaction occurs
            y0 = copy.deepcopy(prev_y0)
            curr_state_after = copy.deepcopy(curr_state)
            curr_state = copy.deepcopy(prev_curr_state)
            floored_curr_state = copy.deepcopy(prev_curr_state)
            propensities_after = copy.deepcopy(propensities)
            propensities = copy.deepcopy(starting_propensities)
            floored_propensities = copy.deepcopy(starting_propensities)
            curr_time = prev_curr_time

            # floored propensites
            for i, s in enumerate(self.model.listOfSpecies):
                floored_curr_state[s] = math.floor(floored_curr_state[s])
            saved_curr_state = copy.deepcopy(curr_state)
            curr_state = floored_curr_state
            for i, r in enumerate(compiled_reactions):
                try:
                    floored_propensities[r] = eval(compiled_reactions[r], {**eval_globals, **curr_state})
                except Exception as e:
                    raise SimulationError('Error calculation propensity for {0}.\nReason: {1}\nfloored_propensities={2}\ncompiled_reactions={3}'.format(r, e, floored_propensities,compiled_reactions))
            curr_state = saved_curr_state

            rxn_times = OrderedDict()
            min_tau = None
            rxn_selected = None
            for rname,rcnt in rxn_count.items():
                if floored_propensities[rname] > 0.0:
                    # estimate the zero crossing time
                    rxn_times[rname] = -1* curr_state[rname] / propensities[rname]
                    if min_tau is None or min_tau > rxn_times[rname]:
                        min_tau = rxn_times[rname]
                        rxn_selected = rname
            if rxn_selected is None: raise SimulationError(f"Negative State detected in step, and no reaction found to fire. error_message={invalid_err_message}")

            tau_step = min_tau #estimated time to the first stochatic reaction

            sol, curr_time = self.__integrate(integrator_options, curr_state,
                                          y0, curr_time, propensities, y_map,
                                          compiled_reactions,
                                          active_rr,
                                          event_queue,
                                          delayed_events,
                                          trigger_states,
                                          event_sensitivity,
                                          tau_step
                                          )

            # only update the selected reaction
            first_rxn_count = copy.deepcopy(rxn_count)
            first_err_message = invalid_err_message
            species_modified,rxn_count = self.__update_stochastic_rxn_states(compiled_reactions, curr_state, only_update=rxn_selected)

            (invalid_state, invalid_err_message) = self.__simulate_invalid_state_check(species_modified, curr_state, compiled_reactions)
            if invalid_state:
                raise SimulationError(f"Negative State detected in step, after single SSA step. error_message={invalid_err_message}")


        save_times, save_index = self.__save_state_to_output(curr_time, save_index, curr_state, species, trajectory, save_times)

        events_processed = self.__process_queued_events(event_queue, trigger_states, curr_state, det_spec)

        # Finally, perform a final check on events after all non-ODE assignment
        # changes have been carried out on model.
        event_cycle = True
        while event_cycle:
            event_cycle = False
            for i, e in enumerate(self.model.listOfEvents.values()):
                triggered = eval(e.trigger.expression, {**eval_globals, **curr_state})
                if triggered and not curr_state[e.name]:
                    curr_state[e.name] = True
                    self.__handle_event(e, curr_state, curr_time,
                                        event_queue, trigger_states, delayed_events)
                    event_cycle = True
                elif not triggered:
                    curr_state[e.name] = False

        events_processed = self.__process_queued_events(event_queue, trigger_states, curr_state, det_spec)

        return sol, curr_state, curr_time, save_times, save_index

    def __set_seed(self, seed):
        # Set seed if supplied
        if seed is not None:
            if not isinstance(seed, int):
                seed = int(seed)
            if seed > 0:
                random.seed(seed)
            else:
                raise ModelError('seed must be a positive integer')

    def __set_recommended_ode_defaults(self, integrator_options):
        """
        Set some ODE solver defaults.  These values are chosen based on the
        precision required to successfully complete the SBML Test suite.
        """

        if 'rtol' not in integrator_options:
            integrator_options['rtol'] = 1e-9
        if 'atol' not in integrator_options:
            integrator_options['atol'] = 1e-12
        if 'max_step' not in integrator_options:
            integrator_options['max_step'] = 0.25

    def __compile_all(self):
        """
        Compile all run-time evaluables to enhance performance.
        """
        compiled_reactions = OrderedDict()
        for i, r in enumerate(self.model.listOfReactions):
            compiled_reactions[r] = compile(self.model.listOfReactions[r].propensity_function, '<string>',
                                            'eval')

        compiled_rate_rules = OrderedDict()
        for i, rr in enumerate(self.model.listOfRateRules.values()):
            if isinstance(rr.variable, str):
                compiled_rate_rules[self.model.listOfSpecies[rr.variable]] = compile(
                                                            rr.formula, '<string>', 'eval')
            else:
                compiled_rate_rules[rr.variable] = compile(rr.formula, '<string>', 'eval')
        compiled_inactive_reactions = OrderedDict()

        compiled_propensities = copy.deepcopy(compiled_reactions)

        return compiled_reactions, compiled_rate_rules, compiled_inactive_reactions, compiled_propensities

    def __initialize_state(self, curr_state, debug):
        """
        Initialize curr_state for each trajectory.
        """

        # intialize parameters to current state
        for p in self.model.listOfParameters:
            curr_state[p] = self.model.listOfParameters[p].value

        # initialize species population state
        for s in self.model.listOfSpecies:
            curr_state[s] = self.model.listOfSpecies[s].initial_value

        # Set reactions to uniform random number
        for i, r in enumerate(self.model.listOfReactions):
            curr_state[r] = math.log(random.uniform(0, 1))
            if debug:
                print("Setting Random number ", curr_state[r], " for ", self.model.listOfReactions[r].name)

        # Initialize event last-fired times to 0
        for e_name in self.model.listOfEvents:
            curr_state[e_name] = 0

        sanitized_species = self.model.sanitized_species_names()
        sanitized_parameters = self.model.sanitized_parameter_names()
        for fd in self.model.listOfFunctionDefinitions.values():
            sanitized_function = fd.sanitized_function(sanitized_species, sanitized_parameters)
            curr_state[fd.name] = eval(f"lambda {', '.join(fd.args)}: {sanitized_function}", eval_globals)

        for ar in self.model.listOfAssignmentRules.values():
            if ar.variable in self.model.listOfSpecies:
                continue
            curr_state[ar.variable] = ar.formula

    def __map_state(self, species, parameters, compiled_reactions, events, curr_state):
        """
        Creates the start state vector for integration and provides a
        dictionary map to it's elements.
        """
        y_map = OrderedDict()
        # Build integration start state
        y0 = [0] * (len(species) + len(parameters) + len(compiled_reactions) + len(events))
        for i, spec in enumerate(species):
            if isinstance(curr_state[spec], str):
                y0[i] = eval(curr_state[spec], {**eval_globals, **curr_state})
            else:
                y0[i] = curr_state[spec]
            y_map[spec] = i
        for i, param in enumerate(parameters):
            y0[i + len(species)] = curr_state[param]
            y_map[param] = i + len(species)
        for i, rxn in enumerate(compiled_reactions):
            y0[i + len(species) + len(parameters)] = curr_state[rxn]
            y_map[rxn] = i + len(species) + len(parameters)
        for i, event in enumerate(events.values()):
            y0[i + len(species) + len(parameters) + len(compiled_reactions)] = curr_state[event.name]
            y_map[event] = i + len(species) + len(parameters) + len(compiled_reactions)
        return y0, y_map

    @classmethod
    def get_solver_settings(cls):
        """
        Returns a list of arguments supported by tau_hybrid_solver.run.
        :returns: Tuple of strings, denoting all keyword argument for this solvers run() method.
        :rtype: tuple
        """
        return ('model', 't', 'number_of_trajectories', 'increment', 'seed', 'debug', 'profile', 'tau_tol',
                'event_sensitivity', 'integrator_options', 'timeout')

    @classmethod
    def get_supported_features(cls):
        return {
            Event,
            RateRule,
            AssignmentRule,
            FunctionDefinition,
        }

    def run(self=None, model=None, t=None, number_of_trajectories=1, increment=None, seed=None,
            debug=False, profile=False, tau_tol=0.03, event_sensitivity=100,integrator_options={},
            live_output=None, live_output_options={}, timeout=None, **kwargs):
        """
        Function calling simulation of the model. This is typically called by the run function in GillesPy2 model
        objects and will inherit those parameters which are passed with the model as the arguments this run function.

        :param model: The model on which the solver will operate. (Deprecated)
        :type model: gillespy2.Model

        :param t: Simulation run time.
        :type t: int or float

        :param number_of_trajectories: Number of trajectories to simulate. By default number_of_trajectories = 1.
        :type number_of_trajectories: int

        :param increment: Save point increment for recording data.
        :type increment: float

        :param seed: The random seed for the simulation. Optional, defaults to None.
        :type seed: int

        :param debug: Set to True to provide additional debug information about the simulation.
        :type debug: bool

        :param profile: Set to True to provide information about step size (tau) taken at each step.
        :type profile: bool

        :param tau_tol: Tolerance level for Tau leaping algorithm.  Larger tolerance values will
        result in larger tau steps. Default value is 0.03.
        :type tau_tol: float

        :param event_sensitivity: Number of data points to be inspected between integration
            steps/save points for event detection. Default event_sensitivity = 100
        :type event_sensitivity: int

        :param integrator_options:  contains options to the scipy integrator. by default, this includes
            rtol=1e-9 and atol=1e-12.  for a list of options,
            see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.LSODA.html.
            Example use: {max_step : 0, rtol : .01}
        :type integrator_options: dict

        :param live_output: The type of output to be displayed by solver. Can be "progress", "text", or "graph".
        :type live_output: str

        :param live_output_options: contains options for live_output. By default {"interval":1}.
            "interval" specifies seconds between displaying.
            "clear_output" specifies if display should be refreshed with each display
        :type live_output_options:  dict

        :param timeout: If set, if simulation takes longer than timeout, will exit.
        :type timeout: int
       
        :returns: A result object containing the results of the simulation.
        :rtype: gillespy2.Results
        """
        from gillespy2 import log

        if self is None:
            # Post deprecation block
            # raise SimulationError("TauHybridSolver must be instantiated to run the simulation")
            # Pre deprecation block
            log.warning(
                """
                `gillespy2.Model.run(solver=TauHybridSolver)` is deprecated.

                You should use `gillespy2.Model.run(solver=TauHybridSolver(model=gillespy2.Model))
                Future releases of GillesPy2 may not support this feature.
                """
            )
            self = TauHybridSolver(model=model)

        if model is not None:
            log.warning('model = gillespy2.model is deprecated. Future releases '
                        'of GillesPy2 may not support this feature.')
        if self.model is None:
            if model is None:
                raise SimulationError("A model is required to run the simulation.")
            self.model = copy.deepcopy(model)

        self.model.compile_prep()
        self.validate_model(self.model, model)
        self.validate_sbml_features(model=self.model)

        self.validate_tspan(increment=increment, t=t)
        if increment is None:
            increment = self.model.tspan[-1] - self.model.tspan[-2]
        if t is None:
            t = self.model.tspan[-1]

        if timeout is not None and timeout > 0:
            for i, s in enumerate(list(self.model._listOfSpecies.keys())):
                # Solve_ivp doesn't return any results until it's finished solving so timing out early only slows
                # the solver.
                if self.model.listOfSpecies[s].mode == 'continuous':
                    timeout = 0
                    log.warning('timeouts not supported by continuous species.')
                    break
                elif self.model.listOfSpecies[s].mode == 'dynamic':
                    log.warning('timeouts not fully supported by dynamic species. If timeout is triggered during'
                                ' integration, total solve time could be longer than expected.')
                    break

        self.stop_event = threading.Event()

        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        if timeout is not None and timeout <= 0:
            timeout = None

        if debug:
            print("t = ", t)
            print("increment = ", increment)

        if len(self.model.listOfEvents):
            self.__set_recommended_ode_defaults(integrator_options)
        self.__set_seed(seed)

        species = list(self.model._listOfSpecies.keys())
        number_species = len(species)

        initial_state = OrderedDict()
        self.__initialize_state(initial_state, debug)
        initial_state['vol'] = self.model.volume
        initial_state['t'] = 0

        # create numpy array for timeline
        timeline = np.linspace(0, t, int(round(t / increment + 1)))
        self.model.tspan = timeline

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline

        # copy initial populations to base
        spec_modes = ['continuous', 'dynamic', 'discrete', None]
        for i, s in enumerate(species):
            if self.model.listOfSpecies[s].mode is None:
                self.model.listOfSpecies[s].mode = 'dynamic'

            if self.model.listOfSpecies[s].mode not in spec_modes:
                raise SpeciesError('Species mode can only be \'continuous\', \'dynamic\',\'discrete\', or '
                                   '\'unspecified(default to dynamic)\'.')
            trajectory_base[:, 0, i + 1] = initial_state[s]

        # curr_time and curr_state are list of len 1 so that __run receives reference
        curr_time = [0]  # Current Simulation Time
        curr_state = [None]
        live_grapher = [None]

        sim_thread = threading.Thread(target=self.___run,
                                      args=(curr_state, curr_time, timeline, trajectory_base, initial_state,
                                            live_grapher,), kwargs={'t': t,
                                                                    'number_of_trajectories': number_of_trajectories,
                                                                    'increment': increment, 'seed': seed,
                                                                    'debug': debug, 'profile': profile,
                                                                    'tau_tol': tau_tol,
                                                                    'event_sensitivity': event_sensitivity,
                                                                    'integrator_options': integrator_options})
        try:
            sim_thread.start()

            if live_output is not None:
                live_output_options['type'] = live_output

                import gillespy2.core.liveGraphing
                gillespy2.core.liveGraphing.valid_graph_params(live_output_options)

                if live_output_options['type'] == "graph":
                    for i, s in enumerate(list(self.model._listOfSpecies.keys())):

                        if self.model.listOfSpecies[s].mode == 'continuous':
                            log.warning('display \"type\" = \"graph\" not recommended with continuous species. '
                                        'Try display \"type\" = \"text\" or \"progress\".')
                            break

                live_grapher[0] = gillespy2.core.liveGraphing.LiveDisplayer(self.model, timeline, number_of_trajectories,
                                                                            live_output_options)
                display_timer = gillespy2.core.liveGraphing.RepeatTimer(live_output_options['interval'],
                                                                        live_grapher[0].display,
                                                                        args=(curr_state, curr_time, trajectory_base, live_output))
                display_timer.start()

            sim_thread.join(timeout=timeout)

            if live_grapher[0] is not None:
                display_timer.cancel()

            self.stop_event.set()
            while self.result is None: pass
        except:
            pass
        if hasattr(self, 'has_raised_exception'):
            raise self.has_raised_exception 

        return Results.build_from_solver_results(self, live_output_options)

    def ___run(self, curr_state, curr_time, timeline, trajectory_base, initial_state, live_grapher, t=20,
               number_of_trajectories=1, increment=0.05, seed=None,
               debug=False, profile=False, tau_tol=0.03, event_sensitivity=100,
               integrator_options={}, **kwargs):
        try:
            self.__run(curr_state, curr_time, timeline, trajectory_base, initial_state, live_grapher, t,
                       number_of_trajectories, increment, seed, debug,
                       profile, tau_tol, event_sensitivity,
                       integrator_options, **kwargs)
        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1

    def __run(self, curr_state, curr_time, timeline, trajectory_base, initial_state, live_grapher, t=20,
              number_of_trajectories=1, increment=0.05, seed=None,
              debug=False, profile=False,
              tau_tol=0.03, event_sensitivity=100,
              integrator_options={}, **kwargs):

        # create mapping of species dictionary to array indices
        species_mappings = self.model._listOfSpecies
        species = list(species_mappings.keys())
        parameter_mappings = self.model._listOfParameters
        parameters = list(parameter_mappings.keys())
        number_species = len(species)

        t0_delayed_events, species_modified_by_events = self.__check_t0_events(initial_state)

        # Create deterministic tracking data structures
        det_spec = {species: value.mode != 'discrete' for (species, value) in self.model.listOfSpecies.items()}

        det_rxn = {rxn: False for (rxn, value) in self.model.listOfReactions.items()}

        # Determine if entire simulation is ODE or Stochastic, in order to
        # avoid unnecessary calculations during simulation
        simulation_data = []

        dependencies = OrderedDict()

        # If considering deterministic changes, create dependency data
        # structure for creating diff eqs later
        for reaction in self.model.listOfReactions:
            dependencies[reaction] = set()
            [dependencies[reaction].add(reactant.name) for reactant in self.model.listOfReactions[reaction].reactants]
            [dependencies[reaction].add(product.name) for product in self.model.listOfReactions[reaction].products]

        # Main trajectory loop
        for trajectory_num in range(number_of_trajectories):

            if self.stop_event.is_set():
                print('exiting')
                self.rc = 33
                break

            # For multi trajectories, live_grapher needs to be informed of trajectory increment
            if live_grapher[0] is not None:
                live_grapher[0].increment_trajectory(trajectory_num)

            trajectory = trajectory_base[trajectory_num]  # NumPy array containing this simulation's results
            propensities = OrderedDict()  # Propensities evaluated at current state

            curr_state[0] = initial_state.copy()
            curr_time[0] = 0  # Current Simulation Time

            end_time = self.model.tspan[-1]  # End of Simulation time
            entry_pos = 1
            data = OrderedDict()  # Dictionary for results
            data['time'] = timeline  # All time entries
            save_times = timeline
            save_index = 0

            # Record Highest Order reactant for each reaction and set error tolerance
            HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, critical_threshold = Tau.initialize(self.model, tau_tol)

            # One-time compilations to reduce time spent with eval
            compiled_reactions, compiled_rate_rules, compiled_inactive_reactions, compiled_propensities = \
                self.__compile_all()
            all_compiled = OrderedDict()
            all_compiled['rxns'] = compiled_reactions
            all_compiled['inactive_rxns'] = compiled_inactive_reactions
            all_compiled['rules'] = compiled_rate_rules

            save_times = np.copy(self.model.tspan)
            delayed_events = []
            trigger_states = {}

            # Handle delayed t0 events
            for state in trigger_states.values():
                if state is None: state = curr_state[0]
            for ename, etime in t0_delayed_events.items():
                curr_state[0][ename] = True
                heapq.heappush(delayed_events, (etime, ename))
                if self.model.listOfEvents[ename].use_values_from_trigger_time:
                    trigger_states[ename] = curr_state[0].copy()
                else:
                    trigger_states[ename] = curr_state[0]

            # Each save step
            while curr_time[0] < self.model.tspan[-1]:

                if self.stop_event.is_set():
                    self.rc = 33
                    break
                # Get current propensities
                for i, r in enumerate(self.model.listOfReactions):
                    try:
                        propensities[r] = eval(compiled_propensities[r], eval_globals, curr_state[0])
                        if curr_state[0][r] > 0 and propensities[r]==0:
                            # This is an edge case, that might happen after a single SSA step.
                            curr_state[0][r] = math.log(random.uniform(0, 1))
                    except Exception as e:
                        raise SimulationError('Error calculation propensity for {0}.\nReason: {1}\ncurr_state={2}'.format(r, e, curr_state))

                # Calculate Tau statistics and select a good tau step
                tau_step = Tau.select(HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, tau_tol, critical_threshold, self.model, propensities, curr_state[0], curr_time[0], save_times[0])

                # Process switching if used
                mn, sd, CV = self.__calculate_statistics(curr_time[0], propensities, curr_state[0], tau_step, det_spec)

                # Calculate sd and CV for hybrid switching and flag deterministic reactions
                deterministic_reactions = self.__flag_det_reactions(det_spec, det_rxn, dependencies)

                if debug:
                    print('mean: {0}'.format(mn))
                    print('standard deviation: {0}'.format(sd))
                    print('CV: {0}'.format(CV))
                    print('det_spec: {0}'.format(det_spec))
                    print('det_rxn: {0}'.format(det_rxn))

                # Set active reactions and rate rules for this integration step
                rr_sets = {frozenset() : compiled_rate_rules} # base rr set
                active_rr = self.__toggle_reactions(all_compiled, deterministic_reactions,
                                                        dependencies, curr_state[0], det_spec, rr_sets)

                # Create integration initial state vector
                y0, y_map = self.__map_state(species, parameters,
                                             compiled_reactions, self.model.listOfEvents, curr_state[0])

                # Run simulation to next step
                sol, curr_state[0], curr_time[0], save_times, save_index = self.__simulate(integrator_options,
                                                                               curr_state[0], y0, curr_time[0],
                                                                               propensities, species,
                                                                               parameters, compiled_reactions,
                                                                               active_rr, y_map,
                                                                               trajectory, save_times, save_index,
                                                                               delayed_events,
                                                                               trigger_states,
                                                                               event_sensitivity, tau_step,
                                                                               debug, det_spec)

            # End of trajectory, format results
            data = {'time': timeline}
            for i in range(number_species):
                data[species[i]] = trajectory[:, i + 1]
            simulation_data.append(data)

        self.result = simulation_data
        return simulation_data, self.rc
