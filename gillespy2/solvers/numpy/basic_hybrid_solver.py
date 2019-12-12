import random, math, sys, warnings
from collections import OrderedDict
from scipy.integrate import ode, solve_ivp
import heapq
import numpy as np
import gillespy2
from gillespy2.core import GillesPySolver, log
from gillespy2.core.gillespyError import *

eval_globals = math.__dict__


class BasicHybridSolver(GillesPySolver):
    """
    This Solver uses a root-finding interpretation of the direct SSA method,
    along with ODE solvers to simulate ODE and Stochastic systems
    interchangeably or simultaneously.
    """
    name = "BasicHybridSolver"

    def __init__(self):
        name = 'BasicHybridSolver'
           
        
    def toggle_reactions(self, model, all_compiled, deterministic_reactions, dependencies, curr_state, det_spec):
        '''
        Helper method which is used to convert reaction channels into
        rate rules, and rate rules into reaction channels, as they are switched
        dynamically throughout the simulation based on user-supplied tolerance.
        '''
        
        #initialize variables
        inactive_reactions = all_compiled['inactive_rxns']
        rate_rules = all_compiled['rules']
        rxns = all_compiled['rxns']
        
        #If the set has changed, reactivate non-determinsitic reactions
        reactivate = []
        for r in inactive_reactions:
            if not r in deterministic_reactions:
                reactivate.append(r)
        for r in reactivate:
            rxns[r] = inactive_reactions.pop(r, None)

        # floor non-det species
        for s, d in det_spec.items():
            if not d and isinstance(curr_state[s], float):
                curr_state[s] = math.floor(curr_state[s])
            
        #Deactivate Determinsitic Reactions
        for r in deterministic_reactions:
            if not r in inactive_reactions:
                inactive_reactions[r] = rxns.pop(r, None)

        #Check if this reaction set is already compiled and in use:
        if deterministic_reactions in rate_rules.keys():
            return

        #Otherwise, this is a new determinstic reaction set that must be compiled
        if not deterministic_reactions in rate_rules:
            rate_rules[deterministic_reactions] = self.create_diff_eqs(deterministic_reactions, model, dependencies)
                
    def create_diff_eqs(self, comb, model, dependencies):
        '''
        Helper method used to convert stochastic reaction descriptions into
        differential equations, used dynamically throught the simulation.
        '''
        diff_eqs = OrderedDict()
        reactions = OrderedDict()
        rate_rules = OrderedDict()

        #Initialize sample dict
        for reaction in comb:
            for dep in dependencies[reaction]:
                if dep not in diff_eqs:
                    diff_eqs[dep] = '0'

        # loop through each det reaction and concatenate it's diff eq for each species
        for reaction in comb:
            factor = OrderedDict()
            for dep in dependencies[reaction]:
                if model.listOfSpecies[dep].mode != 'continuous':
                    pure_continuous = False
            for dep in dependencies[reaction]:
                factor[dep] = 0
            for key, value in model.listOfReactions[reaction].reactants.items():
                factor[key.name] -= value
            for key, value in model.listOfReactions[reaction].products.items():
                factor[key.name] += value
            for dep in dependencies[reaction]:
                if factor[dep] != 0:
                    if model.listOfSpecies[dep].mode == 'continuous':
                        diff_eqs[dep] += ' + {0}*({1})'.format(factor[dep], model.listOfReactions[reaction].ode_propensity_function)
                    else:
                        diff_eqs[dep] += ' + {0}*({1})'.format(factor[dep], model.listOfReactions[reaction].propensity_function)
        
        #create a dictionary of compiled gillespy2 rate rules
        for spec, rate in diff_eqs.items():
            rate_rules[spec] = compile(gillespy2.RateRule(model.listOfSpecies[spec], rate).expression, '<string>', 'eval')

        # be sure to include model rate rules
        for i, rr in enumerate(model.listOfRateRules):
            rate_rules[rr] = compile(model.listOfRateRules[rr].expression, '<string>', 'eval')

        return rate_rules

    def flag_det_reactions(self, model, det_spec, det_rxn, dependencies):
        '''
        Helper method used to flag reactions that can be processed
        deterministically without exceeding the user-supplied tolerance.
        '''
        #Determine if each rxn would be deterministic apart from other reactions
        prev_state = det_rxn.copy()
        for rxn in model.listOfReactions:
            det_rxn[rxn] = True
            for species in dependencies[rxn]:
                if model.listOfSpecies[species].mode == 'discrete': 
                    det_rxn[rxn] = False
                    break
                if model.listOfSpecies[species].mode == 'dynamic' and det_spec[species] == False:
                    det_rxn[rxn] = False
                    break
                    
        # Create a hashable frozenset of all determinstic reactions
        deterministic_reactions = set()
        for rxn in det_rxn:
            if det_rxn[rxn]:
                deterministic_reactions.add(rxn)
        deterministic_reactions = frozenset(deterministic_reactions)
        return deterministic_reactions
                            
    def calculate_statistics(self, *switch_args):
        """
        Calculates Mean, Standard Deviation, and Coefficient of Variance for each
        dynamic species, then set if species can be represented determistically
        """
        mu_i, sigma_i, model, propensities, curr_state, tau_step, det_spec, dependencies, switch_tol = switch_args
        sd = OrderedDict()
        CV = OrderedDict()

        mn = {species:curr_state[species] for (species, value) in 
              model.listOfSpecies.items() if value.mode == 'dynamic'}
        sd = {species:0 for (species, value) in 
              model.listOfSpecies.items() if value.mode == 'dynamic'}

        for r, rxn in model.listOfReactions.items():
                for reactant in rxn.reactants:
                    if reactant.mode == 'dynamic':
                        mn[reactant.name] -= (tau_step * propensities[r] * rxn.reactants[reactant])
                        sd[reactant.name] += (tau_step * propensities[r] * rxn.reactants[reactant]**2)
                for product in rxn.products:
                    if product.mode == 'dynamic':
                        mn[product.name] += (tau_step * propensities[r] * rxn.products[product])
                        sd[product.name] += (tau_step * propensities[r] * rxn.products[product]**2)
                
        # Get coefficient of variance for each dynamic species
        for species in mn:
            if mn[species] > 0:
                CV[species] = sd[species] / mn[species]
            else:
                CV[species] = 1    # value chosen to guarantee discrete
            #Set species to deterministic if CV is less than threshhold
            det_spec[species] = True if CV[species] < switch_tol or model.listOfSpecies[species].mode == 'continuous' else False                            
                
        return sd, CV
    


    @staticmethod
    def __f(t, y, curr_state, species, reactions, rate_rules, propensities,
    y_map, compiled_reactions, compiled_rate_rules, events):
        """
        Evaluate the propensities for the reactions and the RHS of the Reactions and RateRules.
        Also evaluates boolean value of event triggers.
        """
        state_change = [0] * len(y_map)
        curr_state['t'] = t
        for item, index in y_map.items():
            curr_state[item] = y[index]
        for rr in compiled_rate_rules:
            state_change[y_map[rr]] += eval(compiled_rate_rules[rr], eval_globals, curr_state)
        for i, r in enumerate(compiled_reactions):
            propensities[r] = eval(compiled_reactions[r], eval_globals, curr_state)
            state_change[y_map[r]] += propensities[r]
        for event in events:
            triggered = eval(event.trigger.expression, eval_globals, curr_state)
            if triggered: state_change[y_map[event]] = 1


        return state_change

    @staticmethod
    def __event(curr_state, species, reactions, rate_rules, propensities,
    y_map, compiled_reactions, compiled_rate_rules, event_queue, reaction,
    t, y):
        '''
        Base "Event" method used in scipy.integrate.solve_ivp.  This method
        utilizes the brentq method to determine root crossings, and is used in
        conjunction with model stochastic reactions to discover reaction
        firings.
        '''

        return y[y_map[reaction]]
     

    def find_event_time(self, sol, model, start, end, index, depth):
        '''
        Helper method providing binary search implementation for locating
        precise event times.
        '''
        dense_range = np.linspace(start, end, 3)
        mid = dense_range[1]
        if start >= mid or mid >= end or depth == 20: return end
        solutions = np.diff(sol.sol(dense_range)[-len(model.listOfEvents)+index])
        bool_res = [x>0 for x in solutions]

        if bool_res[0]: # event before mid
            depth += 1
            return self.find_event_time(sol, model, dense_range[0],
                dense_range[1], index, depth)
        else: # event after mid
            depth += 1
            return self.find_event_time(sol, model, dense_range[1],
                dense_range[2], index, depth)


    def __detect_events(self, event_sensitivity, sol, model, delayed_events,
                        trigger_states, curr_time):
        '''
        Helper method to locate precise time of event firing.  This method
        first searches for any instance of an event using event_sensitivity to
        determine the granularity of the search.  If an event is detected, a
        binary search is then used to locate the precise time of that event.
        '''
        event_times = {}
        dense_range = np.linspace(sol.t[0], sol.t[-1], len(sol.t)*event_sensitivity)
        for i, e in enumerate(model.listOfEvents.values()):
            solutions = np.diff(sol.sol(dense_range)[-len(model.listOfEvents)+i])
            as_bool = [int(x)>0 for x in solutions]
            bool_res = [x>0 for x in solutions]
            # Search for changes from False to True in event, record first time
            for y in range(len(dense_range)-1):
                # Check Persistent Delays
                if e.name in delayed_events:
                    if not bool_res[y]:
                        heapq.heappop(e.name)
                        del trigger_states[e.name]
                # IF triggered from false to true, refine search
                if bool_res[y] and ((dense_range[y] != curr_time and bool_res[y-1] == 0) or curr_time==0):
                    event_time = self.find_event_time(sol, model, dense_range[y-1],
                        dense_range[y+1], i, 0)
                    if event_time in event_times:
                        event_times[event_time].append(e)
                    else:
                        event_times[event_time] = [e]
                    break
        return event_times
 
    def __get_next_step(self, event_times, reaction_times, delayed_events, sim_end):
        '''
        Helper method to determine the next action to take during simulation,
        and returns that action along with the time that it occurs.
        '''

        next_reaction_time = sim_end + 1
        next_event_trigger = sim_end + 2
        next_delayed_event = sim_end + 3

        curr_time = sim_end
        if len(reaction_times):
            next_reaction_time = min(reaction_times)

        # event triggers
        if (len(event_times)):
            next_event_trigger = min(event_times)

        # delayed events
        if len(delayed_events):
            next_delayed_event = delayed_events[0][0]

        # Execute rxns/events
        next_step = {sim_end: 'end', next_reaction_time: 'rxn', next_event_trigger: 'trigger', next_delayed_event: 'delay'}

        # UPDATE CURR STATE
        curr_time = min(sim_end, next_reaction_time, next_event_trigger,
                        next_delayed_event)
        return next_step[curr_time], curr_time

    def __integrate(self, integrator, integrator_options, curr_state, y0, model, curr_time, 
                                 propensities, y_map, compiled_reactions,
                                 compiled_rate_rules, event_queue,
                                 delayed_events, trigger_states,
                                 event_sensitivity):
        """ 
        Helper function to perform the ODE integration of one step.  This
        method uses scipy.integrate.solve_ivp to get simulation data, and
        determines the next stopping point of the simulation. The state is
        updated and returned to __simulate along with curr_time and the
        solution object. 
        """
        from functools import partial
        events = model.listOfEvents.values()
        int_args = [curr_state, model.listOfSpecies, model.listOfReactions,
                                                          model.listOfRateRules,
                                                          propensities, y_map, 
                                                          compiled_reactions,
                                                          compiled_rate_rules,
                                                          events]

        rhs = lambda t, y: BasicHybridSolver.__f(t, y, *int_args)
        reaction_events = [partial(BasicHybridSolver.__event, *int_args,
        rn) for rn in compiled_reactions]
        for rxn in reaction_events:
            rxn.terminal = True

        curr_state['t'] = curr_time

        # Integrate until end or event is reached
        sol = solve_ivp(rhs, [curr_time, model.tspan[-1]], y0, 
            method=integrator, options=integrator_options, 
            dense_output=True, events=reaction_events, **integrator_options)

        # Search for precise event times
        # TODO MAKE USER INPUT VARIABLE FOR SENSITIVITY
        event_times = self.__detect_events(event_sensitivity, sol, model, delayed_events,
                                    trigger_states, curr_time)

        # rxns
        reaction_times = []
        for rxn in sol.t_events:
            if np.any(rxn):
                reaction_times.append(rxn[0]) # Append first firing time for each rxn

        next_step, curr_time = self.__get_next_step(event_times, reaction_times,
                                                delayed_events, model.tspan[-1])



        for i, s in enumerate(model.listOfSpecies): # Update continuous
            curr_state[s] = sol.sol(curr_time)[i]
        for rxn in compiled_reactions: # Update discrete
            curr_state[rxn] = sol.sol(curr_time)[y_map[rxn]]

        if next_step == 'rxn':
            rxn_state = [sol.sol(curr_time)[y_map[r]] for r in compiled_reactions]
            for rxn in compiled_reactions:
                if curr_state[rxn] == max(rxn_state):
                    curr_state[rxn] = math.log(random.uniform(0, 1))
                    for reactant in model.listOfReactions[rxn].reactants:
                        curr_state[reactant.name] -= model.listOfReactions[rxn].reactants[reactant]
                    for product in model.listOfReactions[rxn].products:
                        curr_state[product.name] += model.listOfReactions[rxn].products[product]
        elif next_step == 'trigger':
            for event in event_times[curr_time]:
                # Fire trigger time events immediately
                if event.delay is None:
                    heapq.heappush(event_queue, (eval(event.priority), event.name))
                # Queue delayed events
                else:
                    curr_state['t'] = curr_time
                    execution_time = curr_time + eval(event.delay,eval_globals, curr_state)
                    heapq.heappush(delayed_events, (execution_time, event.name))
                    trigger_states[event.name] = curr_state.copy()
        elif next_step == 'delay':
            event = heapq.heappop(delayed_events)
            heapq.heappush(event_queue, (eval(model.listOfEvents[event[1]].priority), event[1]))

                

        return sol, curr_time


    def __simulate(self, integrator, integrator_options, curr_state, y0, model, curr_time, 
                        propensities, species, parameters, compiled_reactions,
                        compiled_rate_rules, y_map, trajectory, save_times,
                        delayed_events, trigger_states, event_sensitivity, debug):
        """
        Function to process simulation until next step, which can be a
        stochastic reaction firing, an event trigger or assignment, or end of
        simulation.  Returns three values:
        sol - Python object returned from solve_ivp which contains all solution
        data.
        curr_state - dict containing all state variables for system at current time
        curr_time - float representing current time
        save_times - list of currently unreached save points.
        """

        event_queue = []
        sol, curr_time = self.__integrate(integrator, integrator_options, curr_state, 
                                                           y0, model, curr_time, propensities, y_map, 
                                                           compiled_reactions,
                                                           compiled_rate_rules,
                                                           event_queue,
                                                           delayed_events,
                                                           trigger_states,
                                                           event_sensitivity)


        num_saves = 0
        for time in save_times:
            if time > curr_time: break
            # if a solution is given for it
            trajectory_index = np.where(model.tspan == time)[0][0]
            for s in range(len(species)):
                trajectory[trajectory_index][s+1] = sol.sol(time)[s]
            num_saves += 1
        save_times = save_times[num_saves:]

        pre_assignment_state = curr_state.copy()

        while event_queue:
            # Get events in priority order
            fired_event = model.listOfEvents[heapq.heappop(event_queue)[1]]
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


        return sol, curr_state, curr_time, save_times
    
    def __compile_all(self, model):
        '''
        Compile all run-time evaluables to enhance performance.
        '''
        compiled_reactions = OrderedDict()
        for i, r in enumerate(model.listOfReactions):
            compiled_reactions[r] = compile(model.listOfReactions[r].propensity_function, '<string>',
                                            'eval')
        compiled_rate_rules = OrderedDict()
        for i, rr in enumerate(model.listOfRateRules):
            compiled_rate_rules[rr] = compile(model.listOfRateRules[rr].expression, '<string>', 'eval')
            
        compiled_inactive_reactions = OrderedDict()

        compiled_propensities = OrderedDict()
        for i, r in enumerate(model.listOfReactions):
            compiled_propensities[r] = compile(model.listOfReactions[r].propensity_function, '<string>', 'eval')
        return compiled_reactions, compiled_rate_rules, compiled_inactive_reactions, compiled_propensities

    def __initialize_state(self, model, curr_state, debug):
        '''
        Initialize curr_state for each trajectory.
        '''
        # initialize species population state
        for s in model.listOfSpecies:
            curr_state[s] = model.listOfSpecies[s].initial_value

        # intialize parameters to current state
        for p in model.listOfParameters:
            curr_state[p] = model.listOfParameters[p].value

        # Set reactions to uniform random number
        for i, r in enumerate(model.listOfReactions):
            curr_state[r] = math.log(random.uniform(0, 1))
            if debug:
                print("Setting Random number ", curr_state[r], " for ", model.listOfReactions[r].name)

        # Initialize event last-fired times to 0
        for e_name in model.listOfEvents:
            curr_state[e_name] = 0

    def __map_state(self, species, parameters, compiled_reactions, events, curr_state):
        '''
        Creates the start state vector for integration and provides a
        dictionary map to it's elements.
        '''
        y_map = OrderedDict()
        # Build integration start state
        y0 = [0] * (len(species) + len(parameters) + len(compiled_reactions) + len(events))
        for i, spec in enumerate(species):
            y0[i] = curr_state[spec]
            y_map[spec] = i
        for i, param in enumerate(parameters):
            y0[i+len(species)] = curr_state[param]
            y_map[param] = i+len(species)
        for i, rxn in enumerate(compiled_reactions):
            y0[i+len(species)+len(parameters)] = curr_state[rxn]
            y_map[rxn] = i+len(species)+len(parameters)
        for i, event in enumerate(events.values()):
            y0[i+len(species)+len(parameters)+len(compiled_reactions)] = curr_state[event.name]
            y_map[event] = i+len(species)+len(parameters)+len(compiled_reactions)
        return y0, y_map

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, 
            debug=False, profile=False, show_labels=True, switch_tol=0.03,
            event_sensitivity=100, integrator='LSODA',
            integrator_options={}, **kwargs):
        """
        Function calling simulation of the model. This is typically called by the run function in GillesPy2 model
        objects and will inherit those parameters which are passed with the model as the arguments this run function.

        Attributes
        ----------

        model : GillesPy2.Model
            GillesPy2 model object to simulate
        t : int
            Simulation run time
        number_of_trajectories : int
            The number of times to sample the chemical master equation. Each
            trajectory will be returned at the end of the simulation.
            Optional, defaults to 1.
        increment : float
            Save point increment for recording data
        seed : int
            The random seed for the simulation. Optional, defaults to None.
        debug : bool (False)
            Set to True to provide additional debug information about the
            simulation.
        profile : bool (Fasle)
            Set to True to provide information about step size (tau) taken at each step.
        show_labels: bool (True)
            If true, simulation returns a list of trajectories, where each list entry is a dictionary containing key value pairs of species : trajectory.  If false, returns a numpy array with shape [traj_no, time, species]
        switch_tol: float
            Relative error tolerance value for deterministic/stochastic switching condition between 0.0 and 1.0
        event_sensitivity: int
            Number of data points to be inspected between integration
            steps/save points for event detection
        integrator: String
            integrator method to be used form scipy.integrate.solve_ivp. Options
            include 'RK45', 'RK23', 'Radau', 'BDF', and 'LSODA'.  
            For more details, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
        integrator_options: dictionary
            contains options to the scipy integrator. by default, this includes
            rtol=1e-9 and atol=1e-12.  for a list of options,
            see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html.
            Example use: {max_step : 0, rtol : .01}
        """

        if not isinstance(self, BasicHybridSolver):
            self = BasicHybridSolver()

        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        if debug:
            print("t = ", t)
            print("increment = ", increment)

        if 'rtol' not in integrator_options:
            integrator_options['rtol'] = 1e-9
        if 'atol' not in integrator_options:
            integrator_options['atol'] = 1e-12

        # create mapping of species dictionary to array indices
        species_mappings = model.sanitized_species_names()
        species = list(species_mappings.keys())
        parameter_mappings = model.sanitized_parameter_names()
        parameters = list(parameter_mappings.keys())
        number_species = len(species)

        # create numpy array for timeline
        timeline = np.linspace(0, t, round(t / increment + 1))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline

        spec_modes = ['continuous', 'dynamic', 'discrete']
        # copy initial populations to base
        for i, s in enumerate(species):
            if model.listOfSpecies[s].mode not in spec_modes:
                raise SpeciesError('Species mode can only be \'continuous\', \'dynamic\', or \'discrete\'.')
            trajectory_base[:, 0, i + 1] = model.listOfSpecies[s].initial_value

        det_spec = {species:True for (species, value) in model.listOfSpecies.items() if value.mode == 'dynamic'}
        det_rxn = {rxn:False for (rxn, value) in model.listOfReactions.items()}
        
        dependencies = OrderedDict()

        for reaction in model.listOfReactions:
            dependencies[reaction] = set()
            [dependencies[reaction].add(reactant.name) for reactant in model.listOfReactions[reaction].reactants]
            [dependencies[reaction].add(product.name) for product in model.listOfReactions[reaction].products]

        pure_ode = True
        for reaction in model.listOfReactions.keys():
            for dep in dependencies[reaction]:
                if model.listOfSpecies[dep].mode != 'continuous':
                    pure_ode = False
                    break

        if debug:
            print('dependencies')
            print(dependencies)

        simulation_data = []

        # Set seed if supplied
        if seed is not None:
            if not isinstance(seed, int):
                seed = int(seed)
            if seed > 0:
                random.seed(seed)
            else:
                raise ModelError('seed must be a positive integer')
        for trajectory_num in range(number_of_trajectories):


            trajectory = trajectory_base[trajectory_num] # NumPy array containing this simulation's results
            propensities = OrderedDict() # Propensities evaluated at current state
            curr_state = OrderedDict() # Current state of the system
            curr_time = 0 # Current Simulation Time
            end_time = model.tspan[-1] # End of Simulation time
            entry_pos = 1
            curr_state['vol'] = model.volume # Model volume
            data = OrderedDict() # Dictionary for show_labels results
            data['time'] = timeline # All time entries for show_labels results

            self.__initialize_state(model, curr_state, debug)

            # One-time compilations to reduce time spent with eval
            compiled_reactions, compiled_rate_rules, compiled_inactive_reactions, compiled_propensities = self.__compile_all(model)
            
            all_compiled = OrderedDict()
            all_compiled['rxns'] = compiled_reactions
            all_compiled['inactive_rxns'] = compiled_inactive_reactions
            all_compiled['rules'] = compiled_rate_rules

            save_times = np.copy(model.tspan)
            delayed_events = []
            trigger_states = {}

            # Each save step
            while curr_time < model.tspan[-1]:

                # Get current propensities
                for i, r in enumerate(model.listOfReactions):
                    propensities[r] = eval(compiled_propensities[r], curr_state)

                # Calculate sd and CV for hybrid switching and flag deterministic reactions
                #TODO REWRITE CALCULATION STUFF
                #switch_args = [mu_i, sigma_i, model, propensities, curr_state, tau_step, det_spec, dependencies, switch_tol]
                #sd, CV = self.calculate_statistics(*switch_args)
                deterministic_reactions = self.flag_det_reactions(model, det_spec, det_rxn, dependencies)
                
                if debug:
                    print('mean: {0}'.format(mu_i))
                    print('standard deviation: {0}'.format(sd))
                    print('CV: {0}'.format(CV))
                    print('det_spec: {0}'.format(det_spec))
                    print('det_rxn: {0}'.format(det_rxn))
                
                # Set active reactions and rate rules for this integration step
                self.toggle_reactions(model, all_compiled, deterministic_reactions, dependencies, curr_state, det_spec)
                active_rr = compiled_rate_rules[deterministic_reactions]

                y0, y_map = self.__map_state(species, parameters,
                                        compiled_reactions, model.listOfEvents, curr_state)

                
                #TODO FIX INTEGRATOR OPTIONS

                sol, curr_state, curr_time, save_times = self.__simulate(integrator, integrator_options,
                    curr_state, y0, model, curr_time, propensities, species, 
                    parameters, compiled_reactions, active_rr, y_map,
                    trajectory, save_times, delayed_events, trigger_states,
                    event_sensitivity, debug)

            # End of trajectory
            if show_labels:
                for traj in range(number_of_trajectories):
                    for i in range(number_species):
                        data[species[i]] = trajectory_base[traj, :, i+1]
                    simulation_data.append(data)
            else:
                simulation_data = trajectory_base

            if profile:
                print(steps_taken)
                print("Total Steps Taken: ", len(steps_taken))
                print("Total Steps Rejected: ", steps_rejected)

        return simulation_data
