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

"""Class and methods for the Tau Leaping Solver"""
import copy
import random
import math
from threading import Thread, Event
import numpy as np
from gillespy2.solvers.utilities import Tau
from gillespy2.solvers.utilities import solverutils as nputils
from gillespy2.core import GillesPySolver, log, liveGraphing
from gillespy2.core.gillespyError import *
from gillespy2.core.results import Results

class TauLeapingSolver(GillesPySolver):
    """
    A Tau Leaping solver for GillesPy2 models.  This solver uses an algorithm that calculates
    multiple reactions in a single step over a given tau step size.  The change in propensities
    over this step are bounded by bounding the relative change in state, yielding greatly improved
    run-time performance with very little trade-off in accuracy.

    :param model: The model on which the solver will operate.
    :type model: gillespy2.Model
    """

    name = "TauLeapingSolver"
    rc = 0
    stop_event = None
    pause_event = None
    result = None

    def __init__(self, model=None, debug=False):
        if model is None:
            raise SimulationError("A model is required to run the simulation.")

        name = "TauLeapingSolver"
        rc = 0
        stop_event = None
        pause_event = None
        result = None
        self.model = copy.deepcopy(model)
        self.debug = debug
        self.is_instantiated = True

    def __get_reactions(self, step, curr_state, curr_time, save_time, propensities, reactions):
        """
        Helper Function to get reactions fired from t to t+tau.  
        :returns: Three values:
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
    def get_solver_settings(cls):
        """
        Returns a list of arguments supported by tau_leaping_solver.run.
        :returns: Tuple of strings, denoting all keyword argument for this solvers run() method.
        :rtype: tuple
        """
        return ('model', 't', 'number_of_trajectories', 'increment', 'seed', 'debug', 'profile','timeout', 'tau_tol')

    def run(self=None, model=None, t=None, number_of_trajectories=1, increment=None, seed=None,
            debug=False, profile=False,  live_output=None, live_output_options={},
            timeout=None, resume=None, tau_tol=0.03, **kwargs):
        """
        Function calling simulation of the model.
        This is typically called by the run function in GillesPy2 model objects
        and will inherit those parameters which are passed with the model
        as the arguments this run function.

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

        :param live_output: The type of output to be displayed by solver. Can be "progress", "text", or "graph".
        :type live_output: str

        :param live_output_options: Contains options for live_output. By default {"interval":1}. "interval"
            specifies seconds between displaying. "clear_output" specifies if display should be refreshed with each
            display.
        :type live_output_options: dict

        :param timeout: If set, if simulation takes longer than timeout, will exit.
        :type timeout: int

        :param resume: Result of a previously run simulation, to be resumed.
        :type resume: gillespy2.Results

        :param tau_tol: Tolerance level for Tau leaping algorithm.  Larger tolerance values will
            result in larger tau steps. Default value is 0.03.
        :type tau_tol: float
        
        :returns: A result object containing the results of the simulation.
        :rtype: gillespy2.Results
        """
        from gillespy2 import log

        if self is None:
            # Post deprecation block
            # raise SimulationError("TauLeapingSolver must be instantiated to run the simulation")
            # Pre deprecation block
            log.warning(
                """
                `gillespy2.Model.run(solver=TauLeapingSolver)` is deprecated.

                You should use `gillespy2.Model.run(solver=TauLeapingSolver(model=gillespy2.Model))
                Future releases of GillesPy2 may not support this feature.
                """
            )
            self = TauLeapingSolver(model=model, debug=debug, profile=profile)

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

        self.stop_event = Event()
        self.pause_event = Event()

        if timeout is not None and timeout <= 0:
            timeout = None
        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        # create numpy array for timeline
        if resume is not None:
            # start where we last left off if resuming a simulatio
            lastT = resume['time'][-1]
            step = lastT - resume['time'][-2]
            timeline = np.arange(lastT, t+step, step)
        else:
            timeline = np.linspace(0, t, int(round(t / increment + 1)))

        species = list(self.model._listOfSpecies.keys())
        trajectory_base, tmpSpecies = nputils.numpy_trajectory_base_initialization(self.model, number_of_trajectories,
                                                                                   timeline, species, resume=resume)

        # total_time and curr_state are list of len 1 so that __run receives reference
        if resume is not None:
            total_time = [resume['time'][-1]]
        else:
            total_time = [0]

        curr_state = [None]
        live_grapher = [None]

        sim_thread = Thread(target=self.___run, args=(curr_state, total_time, timeline, trajectory_base, tmpSpecies,
                                                      live_grapher,), kwargs={'t': t,
                                                                              'number_of_trajectories':
                                                                                  number_of_trajectories,
                                                                              'increment': increment, 'seed': seed,
                                                                              'debug': debug, 'resume': resume,
                                                                              'timeout': timeout, 'tau_tol': tau_tol
                                                                              })
        try:
            time = 0
            sim_thread.start()
            if live_output is not None:
                import gillespy2.core.liveGraphing
                live_output_options['type'] = live_output
                gillespy2.core.liveGraphing.valid_graph_params(
                    live_output_options)
                if resume is not None:
                    resumeTest = True  # If resuming, relay this information to live_grapher
                else:
                    resumeTest = False
                live_grapher[
                    0] = gillespy2.core.liveGraphing.LiveDisplayer(self.model,
                                                                   timeline,
                                                                   number_of_trajectories,
                                                                   live_output_options,
                                                                   resume=resumeTest)
                display_timer = gillespy2.core.liveGraphing.RepeatTimer(
                    live_output_options['interval'],
                    live_grapher[0].display, args=(curr_state,
                                                   total_time,
                                                   trajectory_base,
                                                   live_output
                                                   )
                    )
                display_timer.start()

            if timeout is not None:
                while sim_thread.is_alive():
                    sim_thread.join(.1)
                    time += .1
                    if time >= timeout:
                        break
            else:
                while sim_thread.is_alive():
                    sim_thread.join(.1)

            if live_grapher[0] is not None:
                display_timer.cancel()
            self.stop_event.set()
            while self.result is None:
                pass
        except KeyboardInterrupt:
            if live_output:
                display_timer.pause = True
                display_timer.cancel()
            self.pause_event.set()
            while self.result is None:
                pass
        if hasattr(self, 'has_raised_exception'):
            raise SimulationError(
                f"Error encountered while running simulation:\nReturn code: {int(self.rc)}.\n"
            ) from self.has_raised_exception

        return Results.build_from_solver_results(self, live_output_options)

    def ___run(self, curr_state,total_time, timeline, trajectory_base, tmpSpecies, live_grapher, t=20,
               number_of_trajectories=1, increment=0.05, seed=None, debug=False, profile=False,
               resume=None, tau_tol=0.03, **kwargs):

        try:
            self.__run(curr_state, total_time, timeline, trajectory_base, tmpSpecies, live_grapher, t, number_of_trajectories,
                       increment, seed, debug, profile, resume, tau_tol, **kwargs)

        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1

    def __run(self, curr_state, total_time, timeline, trajectory_base, tmpSpecies, live_grapher, t=20,
              number_of_trajectories=1, increment=0.05, seed=None, debug=False, profile=False,
              resume=None, tau_tol=0.03, **kwargs):

        # for use with resume, determines how much excess data to cut off due to
        # how species and time are initialized to 0
        timeStopped = 0
        if resume is not None:
            if resume[0].model != self.model:
                raise ModelError('When resuming, one must not alter the model being resumed.')
            if t < resume['time'][-1]:
                raise ExecutionError(
                    "'t' must be greater than previous simulations end time, or set in the run() method as the "
                    "simulations next end time")
        if debug:
            print("t = ", t)
            print("increment = ", increment)

        species_mappings, species, parameter_mappings, number_species = nputils.numpy_initialization(self.model)

        if seed is not None:
            if not isinstance(seed, int):
                seed = int(seed)
            if seed > 0:
                random.seed(seed)
                np.random.seed(seed)
            else:
                raise ModelError('seed must be a positive integer')

        simulation_data = []

        for trajectory_num in range(number_of_trajectories):
            if self.stop_event.is_set():
                self.rc = 33
                break
            elif self.pause_event.is_set():
                timeStopped = timeline[entry_count]
            # For multi trajectories, live_grapher needs to be informed of trajectory increment
            if live_grapher[0] is not None:
                live_grapher[0].increment_trajectory(trajectory_num)

            start_state = [0] * (len(self.model.listOfReactions) + len(self.model.listOfRateRules))
            propensities = {}
            curr_state[0] = {}

            if resume is not None:
                save_time = total_time[0]
                curr_time = [save_time]
            else:
                save_time = 0
                curr_time = [0]

            curr_state[0]['vol'] = self.model.volume
            data = {'time': timeline}
            steps_taken = []
            steps_rejected = 0
            entry_count = 0
            trajectory = trajectory_base[trajectory_num]

            HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, critical_threshold = Tau.initialize(self.model, tau_tol)

            if resume is not None:
                for spec in self.model.listOfSpecies:
                    curr_state[0][spec] = tmpSpecies[spec]
            else:
                for spec in self.model.listOfSpecies:
                    curr_state[0][spec] = self.model.listOfSpecies[spec].initial_value

            for param in self.model.listOfParameters:
                curr_state[0][param] = self.model.listOfParameters[param].value

            for i, rxn in enumerate(self.model.listOfReactions):
                # set reactions to uniform random number and add to start_state
                start_state[i] = (math.log(random.uniform(0, 1)))
                if debug:
                    print("Setting Random number ",
                          start_state[i], " for ", self.model.listOfReactions[rxn].name)

            compiled_propensities = {}
            for i, r in enumerate(self.model.listOfReactions):
                compiled_propensities[r] = compile(self.model.listOfReactions[r].propensity_function, '<string>', 'eval')

            timestep = 0
            
            # Each save step
            while entry_count < timeline.size:
                if self.stop_event.is_set():
                    self.rc = 33
                    break
                elif self.pause_event.is_set():
                    timeStopped = timeline[entry_count]
                    break

                # Until save step reached
                while curr_time[0] < save_time:
                    if self.stop_event.is_set():
                        self.rc = 33
                        break
                    elif self.pause_event.is_set():
                        timeStopped = timeline[entry_count]
                        break

                    propensity_sum = 0

                    for i, r in enumerate(self.model.listOfReactions):
                        propensities[r] = eval(compiled_propensities[r], curr_state[0])
                        propensity_sum += propensities[r]

                    tau_args = [HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, tau_tol, critical_threshold,
                                self.model, propensities, curr_state[0], curr_time[0], save_time]

                    tau_step = Tau.select(*tau_args)

                    prev_start_state = start_state.copy()
                    prev_curr_state = curr_state[0].copy()
                    prev_curr_time = curr_time[0]
                    total_time[0] = curr_time[0]

                    loop_cnt = 0
                    while True:
                        loop_cnt += 1
                        if loop_cnt > 100:
                            raise Exception("Loop over __get_reactions() exceeded loop count")

                        reactions, curr_state[0], curr_time[0] = self.__get_reactions(
                            tau_step, curr_state[0], curr_time[0], save_time,
                            propensities, self.model.listOfReactions)

                        # Update curr_state with the result of the SSA reaction that fired
                        species_modified = {}
                        for i, rxn in enumerate(self.model.listOfReactions):
                            if reactions[rxn] > 0:
                                for reactant in self.model.listOfReactions[rxn].reactants:
                                    species_modified[reactant.name] = True
                                    curr_state[0][reactant.name] -= self.model.listOfReactions[
                                        rxn].reactants[reactant] * reactions[rxn]
                                for product in self.model.listOfReactions[rxn].products:
                                    species_modified[product.name] = True
                                    curr_state[0][product.name] += self.model.listOfReactions[
                                        rxn].products[product] * reactions[rxn]
                        neg_state = False
                        for spec in species_modified:
                            if curr_state[0][spec] < 0:
                                neg_state = True
                                if debug:
                                    print("Negative state detected: curr_state[{0}]= {1}".format(spec,
                                                                                                 curr_state[0][spec]))
                        if neg_state:
                            if debug:
                                print("\trxn={0}".format(reactions))
                            start_state = prev_start_state.copy()
                            curr_state[0] = prev_curr_state.copy()
                            curr_time[0] = prev_curr_time
                            total_time[0] = prev_curr_time

                            tau_step = tau_step / 2
                            steps_rejected += 1
                            if debug:
                                print("Resetting curr_state[{0}]= {1}".format(spec, curr_state[0][spec]))
                            if debug:
                                print("\tRejecting step, taking step of half size, tau_step={0}".format(tau_step))
                        else:
                            break  # breakout of the while True

                # save step reached
                for i in range(number_species):
                    trajectory[entry_count][i + 1] = curr_state[0][species[i]]
                save_time += increment
                timestep += 1
                entry_count += 1
                
            # end of trajectory
            for i in range(number_species):
                data[species[i]] = trajectory[:, i+1]
            simulation_data.append(data)

            if profile:
                print(steps_taken)
                print("Total Steps Taken: ", len(steps_taken))
                print("Total Steps Rejected: ", steps_rejected)

        # If simulation has been paused, or tstopped !=0
        if timeStopped != 0 or resume is not None:
            simulation_data = nputils.numpy_resume(timeStopped, simulation_data, resume=resume)

        self.result = simulation_data
        return self.result, self.rc
