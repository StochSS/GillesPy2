
import random, math, sys, warnings
from collections import OrderedDict
from scipy.integrate import ode
import heapq
import numpy as np
import threading
import gillespy2
from gillespy2.solvers.numpy import Tau
from gillespy2.core import GillesPySolver, log
from gillespy2.core.gillespyError import *

eval_globals = math.__dict__

def __piecewise(*args):
    '''
    Eval entry for piecewise functions
    '''
    args = list(args)
    sol = None
    if len(args) % 2: args.append(True)
    for i, arg in enumerate(args):
        if not i % 2: continue
        if arg:
            sol = args[i-1]
            break
    return sol
def __xor(*args):
    '''
    Eval entry for MathML xor function
    '''
    from operator import ixor
    from functools import reduce
    args = list(args)
    return reduce(ixor, args)

eval_globals['false'] = False
eval_globals['true'] = True
eval_globals['piecewise'] = __piecewise
eval_globals['xor'] = __xor


class BasicTauHybridSolver(GillesPySolver):
    """
    This Solver uses a root-finding interpretation of the direct SSA method,
    along with ODE solvers to simulate ODE and Stochastic systems
    interchangeably or simultaneously.
    """
    name = "BasicTauHybridSolver"
    rc = 0
    result = None
    stop_event = None

    def __init__(self):
        name = 'BasicTauHybridSolver'
        rc = 0  
        
    def __toggle_reactions(self, model, all_compiled, deterministic_reactions, dependencies, curr_state, det_spec):
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
            rate_rules[deterministic_reactions] = self.__create_diff_eqs(deterministic_reactions, model, dependencies)
                
    def __create_diff_eqs(self, comb, model, dependencies):

        '''
        Helper method used to convert stochastic reaction descriptions into
        differential equations, used dynamically throught the simulation.
        '''
        diff_eqs = OrderedDict()
        rate_rules = OrderedDict()

        #Initialize sample dict
        for spec in model.listOfSpecies:
            if spec in model.listOfRateRules:
                diff_eqs[spec] = model.listOfRateRules[spec].formula
            else:
                diff_eqs[spec] = '0'

        # loop through each det reaction and concatenate it's diff eq for each species
        for reaction in comb:
            factor = {dep:0 for dep in dependencies[reaction]}

            for key, value in model.listOfReactions[reaction].reactants.items():
                if not key.constant and not key.boundary_condition:
                    factor[key.name] -= value
            for key, value in model.listOfReactions[reaction].products.items():
                if not key.constant and not key.boundary_condition:
                    factor[key.name] += value

            for dep in dependencies[reaction]:
                if factor[dep] != 0:
                    if model.listOfSpecies[dep].mode == 'continuous':
                        diff_eqs[dep] += ' + {0}*({1})'.format(factor[dep], model.listOfReactions[reaction].ode_propensity_function)
                    else:
                        diff_eqs[dep] += ' + {0}*({1})'.format(factor[dep], model.listOfReactions[reaction].propensity_function)

        for spec in model.listOfSpecies:
            if diff_eqs[spec] == '0':
                del diff_eqs[spec]
        
        #create a dictionary of compiled gillespy2 rate rules
        for spec, rate in diff_eqs.items():
            rate_rules[spec] = compile(gillespy2.RateRule(model.listOfSpecies[spec], rate).formula, '<string>', 'eval')

        return rate_rules

    def __flag_det_reactions(self, model, det_spec, det_rxn, dependencies):
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
                            
    def __calculate_statistics(self, *switch_args):
        """
        Calculates Mean, Standard Deviation, and Coefficient of Variance for each
        dynamic species, then set if species can be represented determistically
        """
        model, propensities, curr_state, tau_step, det_spec = switch_args

        CV = OrderedDict()
        mn = {species:curr_state[species] for species, value in 
              model.listOfSpecies.items() if value.mode == 'dynamic'}
        sd = {species:0 for species, value in 
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
            sref = model.listOfSpecies[species]
            if sref.switch_min==0:
                if mn[species] > 0:
                    CV[species] = sd[species] / mn[species]
                else:
                    CV[species] = 1    # value chosen to guarantee discrete
                #Set species to deterministic if CV is less than threshhold
                det_spec[species] = CV[species] < sref.switch_tol
            else:
                det_spec[species] = mn[species] > sref.switch_min
        return sd, CV


    @staticmethod
    def __f(t, y, curr_state, species, reactions, rate_rules, propensities,
    y_map, compiled_reactions, compiled_rate_rules, events, assignment_rules):
        """
        Evaluate the propensities for the reactions and the RHS of the Reactions and RateRules.
        Also evaluates boolean value of event triggers.
        """
        state_change = [0] * len(y_map)
        curr_state['t'] = t
        curr_state['time'] = t
        for item, index in y_map.items():
            if item in assignment_rules:
                curr_state[item] = eval(assignment_rules[item].formula, {**eval_globals, **curr_state})
            else:
                curr_state[item] = y[index]
        for rr in compiled_rate_rules:
            try:
                state_change[y_map[rr]] += eval(compiled_rate_rules[rr], {**eval_globals, **curr_state})
            except ValueError:
                pass
        for i, r in enumerate(compiled_reactions):
            propensities[r] = eval(compiled_reactions[r],{**eval_globals, **curr_state})
            state_change[y_map[r]] += propensities[r]
        for event in events:
            triggered = eval(event.trigger.expression, {**eval_globals, **curr_state})
            if triggered: state_change[y_map[event]] = 1


        return state_change


    def __update_stochastic_rxn_states(self, model, compiled_reactions, curr_state):
        '''
        Helper method for updating the state of stochastic reactions.
        '''
        rxn_count = OrderedDict()
        species_modified = OrderedDict()
        # Update stochastic reactions
        for rxn in compiled_reactions:
            rxn_count[rxn] = 0
            while curr_state[rxn] > 0:
                rxn_count[rxn] += 1
                curr_state[rxn] += math.log(random.uniform(0,1))
            if rxn_count[rxn]:
                for reactant in model.listOfReactions[rxn].reactants:
                    species_modified[reactant.name] = True
                    curr_state[reactant.name] -= model.listOfReactions[rxn].reactants[reactant] * rxn_count[rxn]
                for product in model.listOfReactions[rxn].products:
                    species_modified[product.name] = True
                    curr_state[product.name] += model.listOfReactions[rxn].products[product] * rxn_count[rxn]
        return species_modified

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
        '''Set some ODE solver defaults.  These values are chosen based on the
        precision required to successfully complete the SBML Test suite. '''

        if 'rtol' not in integrator_options:
            integrator_options['rtol'] = 1e-9
        if 'atol' not in integrator_options:
            integrator_options['atol'] = 1e-12
        if 'max_step' not in integrator_options:
            integrator_options['max_step'] = 0.25

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
            compiled_rate_rules[rr] = compile(model.listOfRateRules[rr].formula, '<string>', 'eval')
            
        compiled_inactive_reactions = OrderedDict()

        compiled_propensities = compiled_reactions.copy()

        return compiled_reactions, compiled_rate_rules, compiled_inactive_reactions, compiled_propensities

    def __initialize_state(self, model, curr_state, debug):
        '''
        Initialize curr_state for each trajectory.
        '''

        # intialize parameters to current state
        for p in model.listOfParameters:
            curr_state[p] = model.listOfParameters[p].value

        # initialize species population state
        for s in model.listOfSpecies:
            curr_state[s] = model.listOfSpecies[s].initial_value

        # Set reactions to uniform random number
        for i, r in enumerate(model.listOfReactions):
            curr_state[r] = math.log(random.uniform(0, 1))
            if debug:
                print("Setting Random number ", curr_state[r], " for ", model.listOfReactions[r].name)

        # Initialize event last-fired times to 0
        for e_name in model.listOfEvents:
            curr_state[e_name] = 0

        for fd in model.listOfFunctionDefinitions.values():
            curr_state[fd.name] = fd.function

        for ar in model.listOfAssignmentRules.values():
            if ar.variable in model.listOfSpecies: continue
            curr_state[ar.variable] = ar.formula

    def __map_state(self, species, parameters, compiled_reactions, events, curr_state):
        '''
        Creates the start state vector for integration and provides a
        dictionary map to it's elements.
        '''
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
            debug=False, profile=False, show_labels=True,
            tau_tol=0.03, event_sensitivity=100, integrator='LSODA',
            integrator_options={}, timeout=None, **kwargs):
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
        tau_tol: float
            Tolerance level for Tau leaping algorithm.  Larger tolerance values will
            result in larger tau steps. Default value is 0.03.
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

        if isinstance(self, type):
            self = BasicTauHybridSolver()

        self.stop_event = threading.Event()

        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        if timeout is not None and timeout <= 0: timeout = None

        sim_thread = threading.Thread(target=self.___run, args=(model,), kwargs={'t':t,
                                        'number_of_trajectories':number_of_trajectories,
                                        'increment':increment, 'seed':seed,
                                        'debug':debug, 'profile':profile,'show_labels':show_labels,
                                        'timeout':timeout, 'tau_tol':tau_tol,
                                        'event_sensitivity':event_sensitivity,
                                        'integrator':integrator,
                                        'integrator_options':integrator_options})
        try:
            sim_thread.start()
            sim_thread.join(timeout=timeout)
            self.stop_event.set()
            while self.result is None: pass
        except:
            pass
        if hasattr(self,'has_raised_exception'):
            raise self.has_raised_exception
        return self.result, self.rc

    def ___run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, 
            debug=False, profile=False, show_labels=True,
            tau_tol=0.03, event_sensitivity=100, integrator='LSODA',
            integrator_options={}, **kwargs):
            try:
                self.__run(model,t,number_of_trajectories, increment, seed, debug,
                           profile,show_labels, tau_tol, event_sensitivity, integrator,
                           integrator_options, **kwargs)
            except Exception as e:
                self.has_raised_exception = e
                self.result = []
                return [], -1
                
    def __run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, 
            debug=False, profile=False, show_labels=True,
            tau_tol=0.03, event_sensitivity=100, integrator='LSODA',
            integrator_options={}, **kwargs):

        if debug:
            print("t = ", t)
            print("increment = ", increment)

        if len(model.listOfEvents):
            self.__set_recommended_ode_defaults(integrator_options)
        self.__set_seed(seed)

        # create mapping of species dictionary to array indices
        species_mappings = model._listOfSpecies
        species = list(species_mappings.keys())
        parameter_mappings = model._listOfParameters
        parameters = list(parameter_mappings.keys())
        number_species = len(species)

        initial_state = OrderedDict()
        self.__initialize_state(model, initial_state, debug)
        initial_state['vol'] = model.volume
        initial_state['t'] = 0

        # create numpy array for timeline
        timeline = np.linspace(0, t, int(round(t / increment + 1)))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline

        #t0_delayed_events, species_modified_by_events = self.__check_t0_events(model, initial_state)

        # copy initial populations to base
        spec_modes = ['continuous', 'dynamic', 'discrete']
        for i, s in enumerate(species):
            if model.listOfSpecies[s].mode not in spec_modes:
                raise SpeciesError('Species mode can only be \'continuous\', \'dynamic\', or \'discrete\'.')
            trajectory_base[:, 0, i+1] = initial_state[s]

        # Create deterministic tracking data structures
        det_spec = {species:True for (species, value) in model.listOfSpecies.items() if value.mode == 'dynamic'}
        det_rxn = {rxn:False for (rxn, value) in model.listOfReactions.items()}

        # Determine if entire simulation is ODE or Stochastic, in order to
        # avoid unnecessary calculations during simulation
        pure_ode = True
        pure_stochastic = True
        for spec in model.listOfSpecies.values():
            if spec.mode != 'discrete':
                pure_stochastic = False
            if spec.mode != 'continuous':
                pure_ode = False

        if debug:
            print('dependencies')
            print(dependencies)

        simulation_data = []

        dependencies = OrderedDict()

        # If considering deterministic changes, create dependency data
        # structure for creating diff eqs later
        if not pure_stochastic:
            for reaction in model.listOfReactions:
                dependencies[reaction] = set()
                [dependencies[reaction].add(reactant.name) for reactant in model.listOfReactions[reaction].reactants]
                [dependencies[reaction].add(product.name) for product in model.listOfReactions[reaction].products]


        # Main trajectory loop
        for trajectory_num in range(number_of_trajectories):

            if self.stop_event.is_set():
                print('exiting')
                self.rc = 33
                break

            trajectory = trajectory_base[trajectory_num] # NumPy array containing this simulation's results
            propensities = OrderedDict() # Propensities evaluated at current state
            curr_state = initial_state.copy() # Current state of the system
            curr_time = 0 # Current Simulation Time
            end_time = model.tspan[-1] # End of Simulation time
            entry_pos = 1
            data = OrderedDict() # Dictionary for show_labels results
            data['time'] = timeline # All time entries for show_labels results
            save_times = timeline

            # Record Highest Order reactant for each reaction and set error tolerance
            if not pure_ode:
                HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, critical_threshold = Tau.initialize(model, tau_tol)

            # One-time compilations to reduce time spent with eval
            compiled_reactions, compiled_rate_rules, compiled_inactive_reactions, compiled_propensities = self.__compile_all(model)
            all_compiled = OrderedDict()
            all_compiled['rxns'] = compiled_reactions
            all_compiled['inactive_rxns'] = compiled_inactive_reactions
            all_compiled['rules'] = compiled_rate_rules

            save_times = np.copy(model.tspan)
            delayed_events = []
            trigger_states = {}

            # Create integration initial state vector
            y0, y_map = self.__map_state(species, parameters,
                                    compiled_reactions, model.listOfEvents, curr_state)

            events = model.listOfEvents.values()
            # Initialize Integrator
            active_rr = OrderedDict()
            rhs = ode(BasicTauHybridSolver.__f).set_integrator(integrator, **integrator_options)
            rhs.set_initial_value(y0, curr_time)
            rhs.set_f_params(curr_state, model.listOfSpecies, model.listOfReactions,
                                        model.listOfRateRules, propensities, 
                                        y_map, compiled_reactions, active_rr,
                                        events, model.listOfAssignmentRules)


            # Handle delayed t0 events
            for state in trigger_states.values():
                if state is None: state = curr_state
            # TODO Events
            '''
            for ename, etime in t0_delayed_events.items():
                curr_state[ename] = True
                heapq.heappush(delayed_events, (etime, ename))
                if model.listOfEvents[ename].use_values_from_trigger_time:
                    trigger_states[ename] = curr_state.copy()
                else:
                    trigger_states[ename] = curr_state
            '''
            entry_count = 0
            # Each integration step
            save_stop = increment
            while curr_time < model.tspan[-1]:
                save_step = False

                if self.stop_event.is_set(): 
                    self.rc = 33
                    break
                # Get current propensities
                if not pure_ode:
                    for i, r in enumerate(model.listOfReactions):
                        try:
                            propensities[r] = eval(compiled_propensities[r],{**eval_globals, **curr_state})
                        except Exception as e:
                            raise SimulationError('Error calculation propensity for {0}.\nReason: {1}'.format(r, e))

                # Calculate Tau statistics and select a good tau step
                if not pure_ode:
                    tau_args = [HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, tau_tol, critical_threshold,
                            model, propensities, curr_state, curr_time, save_times[0]]
                tau_step = save_times[-1]-curr_time if pure_ode else Tau.select(*tau_args)
                tau_stop = curr_time + tau_step

                # Process switching if used
                if not pure_stochastic and not pure_ode:
                    switch_args = [model, propensities, curr_state, tau_step, det_spec]
                    sd, CV = self.__calculate_statistics(*switch_args)

                # Calculate sd and CV for hybrid switching and flag deterministic reactions
                if pure_stochastic:
                    deterministic_reactions = frozenset() # Empty if non-det
                else:
                    deterministic_reactions = self.__flag_det_reactions(model, det_spec, det_rxn, dependencies)

                # Set active reactions and rate rules for this integration step
                if pure_stochastic:
                    active_rr = compiled_rate_rules 
                else:
                    self.__toggle_reactions(model, all_compiled, deterministic_reactions, dependencies, curr_state, det_spec)
                    active_rr = compiled_rate_rules[deterministic_reactions]
                    rhs.set_f_params(curr_state, model.listOfSpecies, model.listOfReactions,
                                                model.listOfRateRules, propensities, 
                                                y_map, compiled_reactions, active_rr,
                                                events, model.listOfAssignmentRules)
                        
                while curr_time < tau_stop:
                    next_step = min(save_stop, tau_step+curr_time)
                    if increment <= tau_step:
                        save_step = True
                    
                    if debug:
                        print('mean: {0}'.format(mu_i))
                        print('standard deviation: {0}'.format(sd))
                        print('CV: {0}'.format(CV))
                        print('det_spec: {0}'.format(det_spec))
                        print('det_rxn: {0}'.format(det_rxn))
                    
                    # Run simulation to next step
                    loop_count = 0
                    prev_y0 = y0.copy()
                    prev_curr_state = curr_state.copy()
                    prev_curr_time = curr_time
                    while True:
                        loop_count += 1
                        if loop_count > 100:
                            raise Exception("Loop over __integrate() exceeded loop count")
                        y0 = rhs.integrate(next_step)
                        for i, spec in enumerate(model.listOfSpecies):
                            curr_state[spec] = y0[i]
                        for i, rxn in enumerate(compiled_reactions):
                            curr_state[rxn] = y0[y_map[rxn]]
                        species_modified = self.__update_stochastic_rxn_states(model,
                                                            compiled_reactions, curr_state)

                        # Occasionally, a tau step can result in an overly-aggressive
                        # forward step and cause a species population to fall below 0,
                        # which would result in an erroneous simulation.  If this occurs,
                        # back simulation up one step and attempt forward simulation using
                        # a smaller tau step.
                        neg_state = False
                        for s in species_modified.keys():
                            if curr_state[s] < 0:
                                neg_state = True
                        if neg_state:
                            y0 = prev_y0.copy()
                            curr_state = prev_curr_state.copy()
                            curr_time = prev_curr_time
                            tau_step = tau_step / 2
                        else: break 
                    if next_step == save_stop:
                        entry_count += 1
                        save_stop += increment
                        for i, spec in enumerate(model.listOfSpecies):
                            trajectory[entry_count][i+1] = curr_state[spec]
                    curr_time = next_step

                '''
                sol, curr_state, curr_time, save_times = self.__simulate(integrator, integrator_options,
                    curr_state, y0, model, curr_time, propensities, species, 
                    parameters, compiled_reactions, active_rr, y_map,
                    trajectory, save_times, delayed_events, trigger_states,
                    event_sensitivity, tau_step, pure_ode, debug)
                '''
            # End of trajectory, format results
            if show_labels:
                data = {'time':timeline}
                for i in range(number_species):
                    data[species[i]] = trajectory[:, i+1]
                simulation_data.append(data)
            else:
                simulation_data = trajectory_base

            if profile:
                print(steps_taken)
                print("Total Steps Taken: ", len(steps_taken))
                print("Total Steps Rejected: ", steps_rejected)

        self.result = simulation_data
        return simulation_data, self.rc

