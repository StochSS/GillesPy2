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
from threading import Thread, Event
from gillespy2.core.results import Results
from gillespy2.core import GillesPySolver, log
from gillespy2.core.gillespyError import *
from gillespy2.solvers.utilities import solverutils as nputils
import random
import math
import numpy as np
np.set_printoptions(suppress=True)


class NumPySSASolver(GillesPySolver):

    """
    This solver produces simulations of systems via Stochastic Simulation Algorithm.

    :param model: The model on which the solver will operate.
    :type model: gillespy2.Model
    """

    name = "NumPySSASolver"
    rc = 0
    stop_event = None
    result = None
    pause_event = None

    def __init__(self, model=None):
        if model is None:
            raise SimulationError("A model is required to run the simulation.")

        name = 'NumPySSASolver'
        rc = 0
        stop_event = None
        result = None
        pause_event = None
        self.model = copy.deepcopy(model)
        self.is_instantiated = True

    @classmethod
    def get_solver_settings(cls):
        """
        Returns a list of arguments supported by the ssa_solver.run
        :returns: Tuple of strings, denoting all keyword argument for this solvers run() method.
        :rtype: tuple
        """
        return ('model', 't', 'number_of_trajectories', 'increment', 'seed', 'debug', 'timeout')

    def run(self=None, model=None, t=None, number_of_trajectories=1, increment=None, seed=None, debug=False,
            live_output=None, live_output_options={}, timeout=None, resume=None, **kwargs):

        """
        Run the SSA algorithm. Uses a NumPy array for storing results and for generating the timeline.

        :param model: The model on which the solver will operate. (Deprecated)
        :type model: gillespy2.Model
        
        :param t: The end time of the solver.
        :type t: int or float
        
        :param number_of_trajectories: Number of trajectories to simulate. By default number_of_trajectories = 1.
        :type number_of_trajectories: int
            
        :param increment: The time step of the solution.
        :type increment: float
        
        :param seed: The random seed for the simulation. Defaults to None.
        :type seed: int
        
        :param debug: Set to True to provide additional debug information about the simulation.
        :type debug: bool
        
        :param live_output: The type of output to be displayed by solver. Can be "progress", "text", or "graph".
        :type live_output: str
        
        :param live_output_options: dictionary contains options for live_output. By default {"interval":1}.
            "interval" specifies seconds between displaying.
            "clear_output" specifies if display should be refreshed with each display.
        :type live_output_options:  dict
        
        :param timeout: If set, if simulation takes longer than timeout, will exit.
        :type timeout: int
        
        :param resume: Result of a previously run simulation, to be resumed.
        :type resume: gillespy2.Results

        :returns: A result object containing the results of the simulation.
        :rtype: gillespy2.Results
        """
        from gillespy2 import log

        if self is None:
            # Post deprecation block
            # raise SimulationError("NumPySSASolver must be instantiated to run the simulation")
            # Pre deprecation block
            log.warning(
                """
                `gillespy2.Model.run(solver=NumPySSASolver)` is deprecated.

                You should use `gillespy2.Model.run(solver=NumPySSASolver(model=gillespy2.Model))
                Future releases of GillesPy2 may not support this feature.
                """
            )
            self = NumPySSASolver(model=model)

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
            # start where we last left off if resuming a simulation
            lastT = resume['time'][-1]
            step = lastT - resume['time'][-2]
            timeline = np.arange(lastT, t + step, step)
        else:
            timeline = np.linspace(0, t, int(round(t / increment + 1)))

        species = list(self.model._listOfSpecies.keys())

        trajectory_base, tmpSpecies = nputils.numpy_trajectory_base_initialization(self.model, number_of_trajectories,
                                                                                   timeline, species, resume=resume)

        # curr_time and curr_state are list of len 1 so that __run receives reference
        if resume is not None:
            total_time = [resume['time'][-1]]
        else:
            total_time = [0]

        curr_state = [None]
        live_grapher = [None]

        sim_thread = Thread(target=self.___run, args=(curr_state, total_time, timeline, trajectory_base,
                                                      live_grapher,), kwargs={'t': t, 'number_of_trajectories':
                                                                              number_of_trajectories,
                                                                              'increment': increment,
                                                                              'seed': seed, 'debug': debug,
                                                                              'resume': resume, })
        try:
            time = 0
            sim_thread.start()
            if live_output is not None:
                import gillespy2.core.liveGraphing
                live_output_options['type'] = live_output
                gillespy2.core.liveGraphing.valid_graph_params(live_output_options)
                if resume is not None:
                    resumeTest = True  # If resuming, relay this information to live_grapher
                else:
                    resumeTest = False
                live_grapher[0] = gillespy2.core.liveGraphing.LiveDisplayer(self.model, timeline, number_of_trajectories,
                                                                             live_output_options,resume = resumeTest)
                display_timer = gillespy2.core.liveGraphing.RepeatTimer(live_output_options['interval'],
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

    def ___run(self, curr_state, total_time, timeline, trajectory_base, live_grapher, t=20,
               number_of_trajectories=1, increment=0.05, seed=None, debug=False, resume=None,
               ):

        try:
            self.__run(curr_state, total_time, timeline, trajectory_base, live_grapher, t, number_of_trajectories,
                       increment, seed, debug, resume)
        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1

    def __run(self, curr_state, total_time, timeline, trajectory_base, live_grapher, t=20,
              number_of_trajectories=1, increment=0.05, seed=None, debug=False,
              resume=None,):

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

        random.seed(seed)

        species_mappings, species, parameter_mappings, number_species = nputils.numpy_initialization(self.model)

        # create dictionary of all constant parameters for propensity evaluation
        parameters = {'V': self.model.volume}
        for paramName, param in self.model.listOfParameters.items():
            parameters[parameter_mappings[paramName]] = param.value

        # create mapping of reaction dictionary to array indices
        reactions = list(self.model.listOfReactions.keys())

        # create mapping of reactions, and which reactions depend on their reactants/products
        dependent_rxns = nputils.dependency_grapher(self.model, reactions)
        number_reactions = len(reactions)
        propensity_functions = {}

        # create an array mapping reactions to species modified
        species_changes = np.zeros((number_reactions, number_species))

        # pre-evaluate propensity equations from strings:
        for i, reaction in enumerate(reactions):
            # replace all references to species with array indices
            for j, spec in enumerate(species):
                species_changes[i][j] = self.model.listOfReactions[reaction].products.get(self.model.listOfSpecies[spec], 0) \
                                        - self.model.listOfReactions[reaction].reactants.get(self.model.listOfSpecies[spec], 0)
                if debug:
                    print('species_changes: {0},i={1}, j={2}... {3}'.format(species, i, j, species_changes[i][j]))
            propensity_functions[reaction] = [eval('lambda S:' + self.model.listOfReactions[reaction].
                                                   sanitized_propensity_function(species_mappings, parameter_mappings),
                                                   parameters), i]
        if debug:
            print('propensity_functions', propensity_functions)

        # begin simulating each trajectory
        simulation_data = []
        for trajectory_num in range(number_of_trajectories):
            if self.stop_event.is_set():
                self.rc = 33
                break
            elif self.pause_event.is_set():
                timeStopped = timeline[entry_count]
                break

            # For multi trajectories, live_grapher needs to be informed of trajectory increment
            if live_grapher[0] is not None:
                live_grapher[0].increment_trajectory(trajectory_num)

            # copy initial state data
            trajectory = trajectory_base[trajectory_num]
            entry_count = 1
            curr_state[0] = {}
            # curr_time and curr_state are list of len 1 so that __run receives reference
            if resume is not None:
                curr_time = [resume['time'][-1]]
            else:
                curr_time = [0]

            for spec in self.model.listOfSpecies:
                if resume is not None:
                    curr_state[0][spec] = resume[spec][-1]
                else:
                    curr_state[0][spec] = self.model.listOfSpecies[spec].initial_value

            propensity_sums = np.zeros(number_reactions)
            # calculate initial propensity sums
            while entry_count < timeline.size:
                if self.stop_event.is_set():
                    self.rc = 33
                    break
                elif self.pause_event.is_set():
                    timeStopped = timeline[entry_count]
                    break
                # determine next reaction

                species_states = list(curr_state[0].values())

                for i in range(number_reactions):
                    propensity_sums[i] = propensity_functions[reactions[i]][0](species_states)

                    if debug:
                        print('propensity: ', propensity_sums[i])

                propensity_sum = np.sum(propensity_sums)
                if debug:
                    print('propensity_sum: ', propensity_sum)
                # if no more reactions, quit
                if propensity_sum <= 0:
                    trajectory[entry_count:, 1:] = list(species_states)
                    break

                cumulative_sum = random.uniform(0, propensity_sum)
                rand = random.random()
                curr_time[0] += -math.log(rand) / propensity_sum
                total_time[0] += -math.log(rand) / propensity_sum
                if debug:
                    print('cumulative sum: ', cumulative_sum)
                    print('entry count: ', entry_count)
                    print('timeline.size: ', timeline.size)
                    print('curr_time: ', curr_time[0])
                # determine time passed in this reaction

                while entry_count < timeline.size and timeline[entry_count] <= curr_time[0]:
                    if self.stop_event.is_set():
                        self.rc = 33
                        break
                    elif self.pause_event.is_set():
                        timeStopped = timeline[entry_count]
                        break
                    trajectory[entry_count, 1:] = species_states
                    entry_count += 1

                for potential_reaction in range(number_reactions):
                    cumulative_sum -= propensity_sums[potential_reaction]
                    if debug:
                        print('if <=0, fire: ', cumulative_sum)
                    if cumulative_sum <= 0:

                        for i,spec in enumerate(self.model.listOfSpecies):
                            curr_state[0][spec] += species_changes[potential_reaction][i]

                        reacName = reactions[potential_reaction]
                        if debug:
                            print('current state: ', curr_state[0])
                            print('species_changes: ', species_changes)
                            print('updating: ', potential_reaction)

                        species_states = list(curr_state[0].values())
                        for i in dependent_rxns[reacName]['dependencies']:
                            propensity_sums[propensity_functions[i][1]] = propensity_functions[i][0](species_states)

                            if debug:
                                print('new propensity sum: ', propensity_sums[i])
                        break
            data = {
                'time': timeline
            }
            for i in range(number_species):
                data[species[i]] = trajectory[:, i+1]
            simulation_data.append(data)

        # If simulation has been paused, or tstopped !=0
        if timeStopped != 0 or resume is not None:
            simulation_data = nputils.numpy_resume(timeStopped, simulation_data, resume=resume)

        self.result = simulation_data
        return self.result, self.rc
