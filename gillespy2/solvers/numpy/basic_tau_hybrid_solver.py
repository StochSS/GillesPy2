# TODO
'''

the ODE Solution object returned by dense output can be passed an array of data
points.  We can use this, along with including the boolean results of all
events in the passed state vector in order to obtain a reasonable time for
event firings.  The dense output datapoints will require a granularity, which
can be controlled via some user supplied error tolerance
'''



import random, math, sys, warnings
from collections import OrderedDict
from scipy.integrate import ode, solve_ivp
import heapq
import numpy as np
import gillespy2
from gillespy2.core import GillesPySolver, log
from gillespy2.core.gillespyError import *

eval_globals = math.__dict__


class BasicTauHybridSolver(GillesPySolver):
    """
    This Solver uses an algorithm that combines the Tau-Leaping and Hybrid ODE/Stochastic methods.
    A root-finding integration is considered over all reaction channels and continuous species rate
    rules, allowing both continuous and discrete regimes to be considered without partitioning the
    system.  Multiple reactions are fired in a single tau step, and the relative change in propensities
    is bounded by bounding the relative change in the state of the system, resulting in increased
    run-time performance with little accuracy trade-off.
    """
    name = "BasicTauHybridSolver"

    def __init__(self):
        name = 'BasicTauHybridSolver'
           
        
    def toggle_reactions(self, model, all_compiled, deterministic_reactions, dependencies, curr_state, det_spec):
        
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
    y_map, compiled_reactions, compiled_rate_rules, event_queue, event,
    t, y):

        #Put Reactions Here
        pass
     

    def find_event_time(self, sol, model, start, end, index, depth):
        dense_range = np.linspace(start, end, 3)
        mid = dense_range[1]
        if start >= mid or mid >= end or depth == 20: return end
        #print('search times: ', dense_range)
        # print('sol: ', sol.sol(dense_range))
        solutions = np.diff(sol.sol(dense_range)[-len(model.listOfEvents)+index])
        bool_res = [x>0 for x in solutions] # Maybe compare to trigger value?
        #print(bool_res)
        if bool_res[0]: # event before mid
            #print('before mid')
            depth += 1
            return self.find_event_time(sol, model, dense_range[0],
                dense_range[1], index, depth)
        else: # event after mid
            #print('after mid')
            depth += 1
            return self.find_event_time(sol, model, dense_range[1],
                dense_range[2], index, depth)



    def __integrate(self, integrator, integrator_options, curr_state, y0, model, curr_time, 
                                 propensities, y_map, compiled_reactions,
                                 compiled_rate_rules, event_queue):
        """ Helper function to perform the ODE integration of one step """
        from functools import partial
        events = model.listOfEvents.values()
        int_args = [curr_state, model.listOfSpecies, model.listOfReactions,
                                                          model.listOfRateRules,
                                                          propensities, y_map, 
                                                          compiled_reactions,
                                                          compiled_rate_rules,
                                                          events]

        rhs = lambda t, y: BasicTauHybridSolver.__f(t, y, *int_args)
        event_calls = [partial(BasicTauHybridSolver.__event, *int_args,
        event_queue, e) for e in model.listOfEvents.values()]

        #print('curr time: ', curr_time)
        curr_state['t0'] = curr_time
        curr_state['t'] = curr_time
        event_times = {}

        # Get Truth value of triggers before integration
        for e in model.listOfEvents.values():
            e.trigger.value = eval(e.trigger.expression, eval_globals,
                curr_state)

        # Integrate until end or event is reached
        #print('start time: ', curr_time)
        sol = solve_ivp(rhs, [curr_time, model.tspan[-1]], y0, 
            method='Radau', options=integrator_options, 
            dense_output=True)

        # Search for precise event times
        dense_range = np.linspace(sol.t[0], sol.t[-1], 100*len(sol.t))
        for i, e in enumerate(model.listOfEvents.values()):
            solutions = np.diff(sol.sol(dense_range)[-len(model.listOfEvents)+i])
            as_bool = [int(x)>0 for x in solutions]
            bool_res = [x>0 for x in solutions] # Maybe compare to trigger value?
            # Search for changes from False to True in event, record first time
            for y in range(len(dense_range)-1):
                # IF triggered from false to true, refine search
                if bool_res[y] and dense_range[y] != curr_time and bool_res[y-1] == 0:
                    #print('Event: ', e.name, ' found at ', dense_range[y])
                    event_time = self.find_event_time(sol, model, dense_range[y-1],
                        dense_range[y+1], i, 0)
                    if event_time in event_times:
                        event_times[event_time].append(e)
                    else:
                        event_times[event_time] = [e]
                    break
        #print(event_times)

        if (len(event_times)):
            next_event_time = min(event_times)
            if next_event_time:
                curr_time = next_event_time
            for event in event_times[next_event_time]:
                heapq.heappush(event_queue, (eval(event.priority), event.name))
        else: curr_time = model.tspan[-1]
                


        return sol, curr_time


    def __simulate(self, integrator, integrator_options, curr_state, y0, model, curr_time, 
                        propensities, species, parameters, compiled_reactions,
                        compiled_rate_rules, y_map, trajectory, save_times, debug):
        """
        Function to get reactions fired from t to t+tau.  This function solves for root crossings
        of each reaction channel from over tau step, using poisson random number generation
        to calculate distance to the root.  Returns four values:
        rxn_count - dict with key=Reaction channel value=number of times fired
        current - list containing current displacement of each reaction channel for calculating fired reactions.
        curr_state - dict containing all state variables for system at current time
        curr_time - float representing current time
        """

        event_queue = []
        sol, curr_time = self.__integrate(integrator, integrator_options, curr_state, 
                                                           y0, model, curr_time, propensities, y_map, 
                                                           compiled_reactions,
                                                           compiled_rate_rules,
                                                           event_queue)



        num_saves = 0
        for time in save_times:
            if time > curr_time: break
            # if a solution is given for it
            #print('TOP: TIME = ', time)
            trajectory_index = np.where(model.tspan == time)[0][0]
            for s in range(len(species)):
                trajectory[trajectory_index][s+1] = sol.sol(time)[s]
                #print(sol.sol(time)[s])
            num_saves += 1
        save_times = save_times[num_saves:]
        for i, s in enumerate(species):
            curr_state[s] = sol.sol(curr_time)[i]
        #print('at time: ', curr_time)
        #for spec in model.listOfSpecies:
            #print('State of ', spec, ': ', curr_state[spec])
        while event_queue:
            # Get events in priority order
            fired_event = model.listOfEvents[heapq.heappop(event_queue)[1]]
            for a in fired_event.assignments:
                # Get assignment value
                #print('modifying ', a.variable.name)
                assign_value = eval(a.expression, eval_globals, curr_state)
                # Update state of assignment variable
                curr_state[a.variable.name] = assign_value
        #for spec in model.listOfSpecies:
            #print('State of ', spec, ': ', curr_state[spec])

        # TODO THIS NEEDS TO BE HANDLED IN AN EVENT-LIKE MANNER
        # TODO WILL ALSO NEED TO RECALCULATE PROPENSITIES
        '''
        # UPDATE THE STATE of the discrete reactions
        for i, r in enumerate(compiled_reactions):
            curr_state[r] = current[i+len(species)+len(parameters)]
        rxn_count = OrderedDict()
        fired = False
        for i, r in enumerate(compiled_reactions):
            rxn_count[r] = 0
            while curr_state[r] > 0:
                if not fired:
                    fired = True
                rxn_count[r] += 1
                urn = (math.log(random.uniform(0, 1)))
                current[i+len(compiled_rate_rules)] += urn
                curr_state[r] += urn
        '''

        return sol, curr_state, curr_time, save_times

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, 
            debug=False, profile=False, show_labels=True, switch_tol=0.03, tau_tol=0.03, 
            integrator='RK45', integrator_options={}, **kwargs):
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
        tau_tol: float
            Relative error tolerance value for calculating tau step between 0.0 and 1.0
        integrator: String
            integrator to be used form scipy.integrate.ode. Options include 'vode', 'zvode', 'lsoda', 'dopri5', and 'dop835'.  For more details, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
        integrator_options: dictionary
            contains options to the scipy integrator. for a list of options, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html.
            Example use: {max_step : 0, rtol : .01}
        """

        if not isinstance(self, BasicTauHybridSolver):
            self = BasicTauHybridSolver()

        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        if debug:
            print("t = ", t)
            print("increment = ", increment)

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

            # One-time compilations to reduce time spent with eval
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
            
            all_compiled = OrderedDict()
            all_compiled['rxns'] = compiled_reactions
            all_compiled['rules'] = compiled_rate_rules
            all_compiled['inactive_rxns'] = compiled_inactive_reactions

            save_times = np.copy(model.tspan)

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

                y_map = OrderedDict()
                # Build integration start state
                y0 = [0] * (len(species) + len(parameters) + len(compiled_reactions) + len(model.listOfEvents))
                for i, spec in enumerate(species):
                    y0[i] = curr_state[spec]
                    y_map[spec] = i
                for i, param in enumerate(parameters):
                    y0[i+len(species)] = curr_state[param]
                    y_map[param] = i+len(species)
                for i, rxn in enumerate(compiled_reactions):
                    y0[i+len(species)+len(parameters)] = curr_state[rxn]
                    y_map[rxn] = i+len(species)+len(parameters)
                for i, event in enumerate(model.listOfEvents.values()):
                    y0[i+len(species)+len(parameters)+len(compiled_reactions)] = curr_state[event.name]
                    y_map[event] = i+len(species)+len(parameters)+len(compiled_reactions)

                
                #TODO FIX INTEGRATOR OPTIONS

                sol, curr_state, curr_time, save_times = self.__simulate(integrator, integrator_options,
                    curr_state, y0, model, curr_time, propensities, species, 
                    parameters, compiled_reactions, active_rr, y_map,
                    trajectory, save_times, debug)
                # TODO WILL HAVE TO UPDATE THIS TO NEW API
                '''
                # Update curr_state with the result of the SSA reaction that fired
                species_modified = OrderedDict()
                for i, r in enumerate(compiled_reactions):
                    if reactions[r] > 0:
                        for reactant in model.listOfReactions[r].reactants:
                            species_modified[str(reactant)] = True
                            curr_state[str(reactant)] -= model.listOfReactions[r].reactants[reactant] * reactions[r]
                        for product in model.listOfReactions[r].products:
                            species_modified[str(product)] = True
                            curr_state[str(product)] += model.listOfReactions[r].products[product] * reactions[r]
                '''         


            # TODO FIX SHOW_LABELS
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
