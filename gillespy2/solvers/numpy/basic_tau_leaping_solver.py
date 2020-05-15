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
    pause_event = None
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
        pause_event = None
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
                    debug=False, profile=False, show_labels=True, 
                    timeout=None, resume=None, resumeTime=None, tau_tol=0.03, **kwargs):
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
                resume : Result of previous simulation
                    resultResult of a previously run simulation, to be resumed
                resumeTime : int
                    How much longer to run the previously ran simulation
                """

        if isinstance(self, type):
            self = BasicTauLeapingSolver(debug=debug, profile=profile)

        self.stop_event = Event()
        self.pause_event = Event()

        if timeout is not None and timeout <= 0: timeout = None
        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        sim_thread = Thread(target=self.___run, args=(model,), kwargs={'t':t,
                                        'number_of_trajectories':number_of_trajectories,
                                        'increment':increment, 'seed':seed,
                                        'debug':debug, 'show_labels':show_labels,
                                        'resume':resume, 'resumeTime':resumeTime,
                                        'timeout':timeout, 'tau_tol':tau_tol})
        try:
            sim_thread.start()
            sim_thread.join(timeout=timeout)
            self.stop_event.set()
            while self.result is None: pass
        except KeyboardInterrupt:
            print('interrupted!')
            self.pause_event.set()
            while self.result is None: pass
        if hasattr(self, 'has_raised_exception'):
            raise self.has_raised_exception
        return self.result, self.rc

    def ___run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
                    debug=False, profile=False, show_labels=True, 
                    timeout=None, resume=None, resumeTime=None, tau_tol=0.03, **kwargs):

        if resume != None and resumeTime == None:
            log.warning("If resuming a simulation, must set a 'resumeTime' in the run() function")
        if resumeTime != None and resumeTime<resume['time'][-1]:
            log.warning("resumeTime must be greater than previous simulations end time.")
        elif resume != None and resumeTime != None:
            t = resumeTime

        try:
            self.__run(model, t, number_of_trajectories, increment, seed,
                        debug, profile, show_labels, timeout, resume,resumeTime, tau_tol, **kwargs)
        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1


    def __run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
                    debug=False, profile=False, show_labels=True, 
                    timeout=None, resume=None, resumeTime=None, tau_tol=0.03, **kwargs):

        #For use if pausing a simulation, at what time simulation is stopped
        timeStopped = 0

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

        # create numpy array for timeline
        if resumeTime != None:
            # start where we last left off if resuming a simulation
            timeline = np.linspace(resume['time'][-1], t, int(round(t - resume['time'][-1] + 1)))
        else:
            timeline = np.linspace(0, t, int(round(t / increment + 1)))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline
        # copy initial populations to base

        if resume != None:
            tmpSpecies = {}
            # Set initial values of species to where last left off
            for i in species:
                tmpSpecies[i] = resume[i][-1]
            for i, s in enumerate(species):
                trajectory_base[:, 0, i + 1] = tmpSpecies[s]
        else:
            for i, s in enumerate(species):
                trajectory_base[:, 0, i + 1] = model.listOfSpecies[s].initial_value
                # create dictionary of all constant parameters for propensity evaluation

        simulation_data = []

        for trajectory_num in range(number_of_trajectories):
            if self.stop_event.is_set():
                self.rc = 33
                break
            elif self.pause_event.is_set():
                print(timeStopped)
                timeStopped = timeline[entry_count]

            start_state = [0] * (len(model.listOfReactions) + len(model.listOfRateRules))
            propensities = {}
            curr_state = {}
            curr_time = 0
            curr_state['vol'] = model.volume
            save_time = 0
            data = { 'time': timeline}
            steps_taken = []
            steps_rejected = 0
            entry_count = 0
            trajectory = trajectory_base[trajectory_num]


            HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, critical_threshold = Tau.initialize(model, tau_tol)
            # initialize populations
            if resume != None:
                for spec in model.listOfSpecies:
                    curr_state[spec] = tmpSpecies[spec]
            else:
                for spec in model.listOfSpecies:
                    curr_state[spec] = model.listOfSpecies[spec].initial_value

            for param in model.listOfParameters:
                curr_state[param] = model.listOfParameters[param].value

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
                elif self.pause_event.is_set():
                    timeStopped = timeline[entry_count]
                    print("Time stopped at : "+ str(timeStopped))
                    break
                
                #Until save step reached
                while curr_time < save_time:
                    if self.stop_event.is_set():
                        self.rc = 33
                        break
                    elif self.pause_event.is_set():
                        timeStopped = timeline[entry_count]
                        print("Time stopped at : "+ str(timeStopped))
                        break

                    propensity_sum = 0

                    for i, r in enumerate(model.listOfReactions):
                        propensities[r] = eval(compiled_propensities[r], curr_state)
                        propensity_sum += propensities[r]

                    tau_args = [HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, tau_tol, critical_threshold,
                            model, propensities, curr_state, curr_time, save_time]

                    tau_step = Tau.select(*tau_args)

                    prev_start_state = start_state.copy()
                    prev_curr_state = curr_state.copy()
                    prev_curr_time = curr_time

                    loop_cnt = 0
                    while True:
                        loop_cnt += 1
                        if loop_cnt > 100:
                            raise Exception("Loop over __get_reactions() exceeded loop count")

                        reactions, curr_state, curr_time = self.__get_reactions(
                            tau_step, curr_state, curr_time, save_time,
                            propensities, model.listOfReactions)

                        # Update curr_state with the result of the SSA reaction that fired
                        species_modified = {}
                        for i, rxn in enumerate(model.listOfReactions):
                            if reactions[rxn] > 0:
                                for reactant in model.listOfReactions[rxn].reactants:
                                    species_modified[reactant.name] = True
                                    curr_state[reactant.name] -= model.listOfReactions[
                                        rxn].reactants[reactant] * reactions[rxn]
                                for product in model.listOfReactions[rxn].products:
                                    species_modified[product.name] = True
                                    curr_state[product.name] += model.listOfReactions[
                                        rxn].products[product] * reactions[rxn]
                        neg_state = False
                        for spec in species_modified:
                            if curr_state[spec] < 0:
                                neg_state = True
                                if debug:
                                    print("Negative state detected: curr_state[{0}]= {1}".format(
                                        spec, curr_state[spec]))
                        if neg_state:
                            if debug:
                                print("\trxn={0}".format(reactions))
                            start_state = prev_start_state.copy()
                            curr_state = prev_curr_state.copy()
                            curr_time = prev_curr_time
                            tau_step = tau_step / 2
                            steps_rejected += 1
                            if debug:
                                print("Resetting curr_state[{0}]= {1}".format(
                                    spec, curr_state[spec]))
                            if debug:
                                print(
                                    "\tRejecting step, taking step of half size, ",
                                    "tau_step={0}".format(tau_step))
                        else:
                            break  # breakout of the while True

                # save step reached
                for i in range(number_species):
#                     print('appending {0} to species {1}'.format(curr_state[species[i]], species[i]))
                    trajectory[entry_count][i + 1] = curr_state[species[i]]
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

        # If simulation has been paused, or tstopped !=0

        if timeStopped != 0 and timeStopped != simulation_data[0]['time'][-1]:
            tester = np.where(simulation_data[0]['time'] > timeStopped)[0].size
            index = np.where(simulation_data[0]['time'] == timeStopped)[0][0]
            if tester > 0:
                for i in simulation_data[0]:
                    simulation_data[0][i] = simulation_data[0][i][:index]

        if resume != None:
            # If resuming, combine old pause with new data
            for i in simulation_data[0]:
                oldData = resume[i][:-1]
                newData = simulation_data[0][i]
                simulation_data[0][i] = np.concatenate((oldData, newData), axis=None)

        self.result = simulation_data
        return simulation_data, self.rc
