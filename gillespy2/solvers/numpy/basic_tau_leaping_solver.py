"""Class and methods for Basic Tau Leaping Solver"""


import random, math, sys, warnings
from threading import Thread, Event
import numpy as np
from gillespy2.solvers.numpy import Tau
from gillespy2.core import GillesPySolver, log


class BasicTauLeapingSolver(GillesPySolver):
    name = 'BasicTauLeapingSolver'
    rc = 0
    stop_event = None
    result = None
    """
    A Basic Tau Leaping Solver for GillesPy2 models.  This solver uses an algorithm calculates
    multiple reactions in a single step over a given tau step size.  The change in propensities
    over this step are bounded by bounding the relative change in state, yielding greatly improved
    run-time performance with very little trade-off in accuracy.
    """
    name = "BasicTauLeapingSolver"

    def __init__(self, debug=False, profile=False):
        name = "BasicTauLeapingSolver"
        rc = 0
        stop_event = None
        result = None
        self.debug = debug
        self.profile = profile

    def __get_reactions(self, step, curr_state, curr_time, save_time, propensities, reactions):
        """
        Helper Function to get reactions fired from t to t+tau.  Returns three values:
        rxn_count - dict with key=Reaction channel value=number of times fired
        curr_state - dict containing all state variables for system at current time
        curr_time - float representing current time
        """


        if curr_time + step > save_time:
            if self.debug:
                print("Step exceeds save_time, changing step size from ", step,
                      " to ", save_time - curr_time)
            step = save_time - curr_time

        if self.debug:
            print("Curr Time: ", curr_time, " Save time: ", save_time, "step: ", step)

        rxn_count = {}

        for rxn in reactions:
            rxn_count[rxn] = np.random.poisson(propensities[rxn] * step)

        if self.debug:
            print("Reactions Fired: ", rxn_count)

        curr_time = curr_time+step

        return rxn_count, curr_state, curr_time

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
                    debug=False, profile=False, show_labels=True, display_interval = 0, display_type =None,
                    timeout=None, tau_tol=0.03, **kwargs):
        """
        Function calling simulation of the model.
        This is typically called by the run function in GillesPy2 model objects
        and will inherit those parameters which are passed with the model
        as the arguments this run function.

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
                show_labels : bool (True)
                    Use names of species as index of result object rather than position numbers.
                """

        if isinstance(self, type):
            self = BasicTauLeapingSolver(debug=debug, profile=profile)

        self.stop_event = Event()
        if timeout is not None and timeout <= 0: timeout = None
        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        # create numpy array for timeline
        timeline = np.linspace(0, t, int(round(t / increment + 1)))

        species = list( model._listOfSpecies.keys())
        number_species = len(species)

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline

        # curr_time and curr_state are list of len 1 so that __run receives reference
        curr_time = [0]  # Current Simulation Time
        curr_state = [None]
        live_grapher = [None]

        sim_thread = Thread(target=self.___run, args=(model, curr_state,curr_time, timeline, trajectory_base,live_grapher,), kwargs={'t':t,
                                        'number_of_trajectories':number_of_trajectories,
                                        'increment':increment, 'seed':seed,
                                        'debug':debug, 'show_labels':show_labels,
                                        'timeout':timeout, 'tau_tol':tau_tol})
        try:
            sim_thread.start()

            from gillespy2.core.liveGraphing import valid_graph_params
            if valid_graph_params(display_type, display_interval):
                import gillespy2.core.liveGraphing
                live_grapher[0] = gillespy2.core.liveGraphing.LiveDisplayer(display_type, display_interval, model,
                                                                            timeline.size, number_of_trajectories)
                display_timer = gillespy2.core.liveGraphing.RepeatTimer(display_interval, live_grapher[0].display,
                                                                        args=(curr_state, curr_time, trajectory_base,))
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
        return self.result, self.rc

    def ___run(self, model, curr_state,curr_time, timeline, trajectory_base, live_grapher, t=20, number_of_trajectories=1, increment=0.05, seed=None,
                    debug=False, profile=False, show_labels=True, 
                    timeout=None, tau_tol=0.03, **kwargs):
        try:
            self.__run(model, curr_state,curr_time, timeline, trajectory_base, live_grapher, t, number_of_trajectories, increment, seed,
                        debug, profile, show_labels, timeout, tau_tol, **kwargs)
        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1


    def __run(self, model, curr_state,curr_time, timeline, trajectory_base, live_grapher, t=20, number_of_trajectories=1, increment=0.05, seed=None,
                    debug=False, profile=False, show_labels=True, 
                    timeout=None, tau_tol=0.03, **kwargs):

        if debug:
            print("t = ", t)
            print("increment = ", increment)
            
        species_mappings = model.sanitized_species_names()
        species = list(species_mappings.keys())
        parameter_mappings = model.sanitized_parameter_names()
        number_species = len(species)

        if seed is not None:
            if not isinstance(seed, int):
                seed = int(seed)
            if seed > 0:
                random.seed(seed)
                np.random.seed(seed)
            else:
                raise ModelError('seed must be a positive integer')

        # copy initial populations to base
        for i, s in enumerate(species):
            trajectory_base[:, 0, i + 1] = model.listOfSpecies[s].initial_value
            # create dictionary of all constant parameters for propensity evaluation

        simulation_data = []

        for trajectory_num in range(number_of_trajectories):
            if self.stop_event.is_set():
                self.rc = 33
                break

            # For multi trajectories, live_grapher needs to be informed of trajectory increment
            if live_grapher[0] is not None:
                live_grapher[0].increment_trajectory(trajectory_num)

            start_state = [0] * (len(model.listOfReactions) + len(model.listOfRateRules))
            propensities = {}
            curr_state[0] = {}
            curr_time[0] = 0
            curr_state[0]['vol'] = model.volume
            save_time = 0
            data = { 'time': timeline}
            steps_taken = []
            steps_rejected = 0
            entry_count = 0
            trajectory = trajectory_base[trajectory_num]


            HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, critical_threshold = Tau.initialize(model, tau_tol)

            for spec in model.listOfSpecies:
                # initialize populations
                curr_state[0][spec] = model.listOfSpecies[spec].initial_value

            for param in model.listOfParameters:
                curr_state[0][param] = model.listOfParameters[param].value

            for i, rxn in enumerate(model.listOfReactions):
                # set reactions to uniform random number and add to start_state
                start_state[i] = (math.log(random.uniform(0, 1)))
                if debug:
                    print("Setting Random number ",
                          start_state[i], " for ", model.listOfReactions[rxn].name)

            compiled_propensities = {}
            for i, r in enumerate(model.listOfReactions):
                compiled_propensities[r] = compile(model.listOfReactions[r].propensity_function, '<string>', 'eval')

            timestep = 0
            
            #Each save step
            while entry_count < timeline.size:
                if self.stop_event.is_set():
                    self.rc = 33
                    break
                
                #Until save step reached
                while curr_time[0] < save_time:
                    if self.stop_event.is_set():
                        self.rc = 33
                        break
                    propensity_sum = 0

                    for i, r in enumerate(model.listOfReactions):
                        propensities[r] = eval(compiled_propensities[r], curr_state[0])
                        propensity_sum += propensities[r]

                    tau_args = [HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, tau_tol, critical_threshold,
                            model, propensities, curr_state[0], curr_time[0], save_time]

                    tau_step = Tau.select(*tau_args)

                    prev_start_state = start_state.copy()
                    prev_curr_state = curr_state[0].copy()
                    prev_curr_time = curr_time[0]

                    loop_cnt = 0
                    while True:
                        loop_cnt += 1
                        if loop_cnt > 100:
                            raise Exception("Loop over __get_reactions() exceeded loop count")

                        reactions, curr_state[0], curr_time[0] = self.__get_reactions(
                            tau_step, curr_state[0], curr_time[0], save_time,
                            propensities, model.listOfReactions)

                        # Update curr_state with the result of the SSA reaction that fired
                        species_modified = {}
                        for i, rxn in enumerate(model.listOfReactions):
                            if reactions[rxn] > 0:
                                for reactant in model.listOfReactions[rxn].reactants:
                                    species_modified[reactant.name] = True
                                    curr_state[0][reactant.name] -= model.listOfReactions[
                                        rxn].reactants[reactant] * reactions[rxn]
                                for product in model.listOfReactions[rxn].products:
                                    species_modified[product.name] = True
                                    curr_state[0][product.name] += model.listOfReactions[
                                        rxn].products[product] * reactions[rxn]
                        neg_state = False
                        for spec in species_modified:
                            if curr_state[0][spec] < 0:
                                neg_state = True
                                if debug:
                                    print("Negative state detected: curr_state[{0}]= {1}".format(
                                        spec, curr_state[0][spec]))
                        if neg_state:
                            if debug:
                                print("\trxn={0}".format(reactions))
                            start_state = prev_start_state.copy()
                            curr_state[0] = prev_curr_state.copy()
                            curr_time[0] = prev_curr_time
                            tau_step = tau_step / 2
                            steps_rejected += 1
                            if debug:
                                print("Resetting curr_state[{0}]= {1}".format(
                                    spec, curr_state[0][spec]))
                            if debug:
                                print(
                                    "\tRejecting step, taking step of half size, ",
                                    "tau_step={0}".format(tau_step))
                        else:
                            break  # breakout of the while True

                # save step reached
                for i in range(number_species):
#                     print('appending {0} to species {1}'.format(curr_state[species[i]], species[i]))
                    trajectory[entry_count][i + 1] = curr_state[0][species[i]]
                save_time += increment
                timestep += 1
                entry_count += 1
                
            # end of trajectory
            if show_labels:
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
