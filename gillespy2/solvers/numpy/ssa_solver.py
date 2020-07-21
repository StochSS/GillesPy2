from threading import Thread, Event
from gillespy2.core import GillesPySolver, log, gillespyError
from gillespy2.solvers.utilities import solverutils as nputils
import random
import math
import numpy as np
np.set_printoptions(suppress=True)


class NumPySSASolver(GillesPySolver):
    name = "NumPySSASolver"
    rc = 0
    stop_event = None
    result = None
    pause_event = None

    def __init__(self):
        name = 'NumPySSASolver'
        rc = 0
        stop_event = None
        result = None
        pause_event = None

    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'increment', 'seed', 'debug', 'timeout')

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=True,
            live_output=None, live_output_options={}, timeout=None, resume=None, **kwargs):

        """
        Run the SSA algorithm using a NumPy for storing the data in arrays and generating the timeline.
        :param model: The model on which the solver will operate.
        :param t: The end time of the solver.
        :param number_of_trajectories: The number of times to sample the chemical master equation. Each
        trajectory will be returned at the end of the simulation.
        :param increment: The time step of the solution.
        :param seed: The random seed for the simulation. Defaults to None.
        :param debug: Set to True to provide additional debug information about the
        simulation.
        :param resume: Result of a previously run simulation, to be resumed
        :param live_output : str The type of output to be displayed by solver. Can be "progress", "text", or "graph".
        :param live_output_options : dictionary contains options for live_output. By default {"interval":1}.
                    "interval" specifies seconds between displaying.
                    "clear_output" specifies if display should be refreshed with each display
        :return: a list of each trajectory simulated.
        """

        if isinstance(self, type):
            self = NumPySSASolver()

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

        species = list(model._listOfSpecies.keys())

        trajectory_base, tmpSpecies = nputils.numpy_trajectory_base_initialization(model, number_of_trajectories,
                                                                                   timeline, species, resume=resume)

        # curr_time and curr_state are list of len 1 so that __run receives reference
        if resume is not None:
            total_time = [resume['time'][-1]]
        else:
            total_time = [0]

        curr_state = [None]
        live_grapher = [None]

        sim_thread = Thread(target=self.___run, args=(model, curr_state, total_time, timeline, trajectory_base,
                                                      live_grapher,), kwargs={'t': t, 'number_of_trajectories':
                                                                              number_of_trajectories,
                                                                              'increment': increment,
                                                                              'seed': seed, 'debug': debug,
                                                                              'show_labels': show_labels,
                                                                              'timeout': timeout,
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
                live_grapher[0] = gillespy2.core.liveGraphing.LiveDisplayer(model, timeline, number_of_trajectories,
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
            raise self.has_raised_exception

        return self.result, self.rc

    def ___run(self, model, curr_state, total_time, timeline, trajectory_base, live_grapher, t=20,
               number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=True, resume=None,
               timeout=None):

        try:
            self.__run(model, curr_state, total_time, timeline, trajectory_base, live_grapher, t, number_of_trajectories,
                       increment, seed, debug, show_labels, resume, timeout)
        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1

    def __run(self, model, curr_state, total_time, timeline, trajectory_base, live_grapher, t=20,
              number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=True,
              resume=None,  timeout=None):

        # for use with resume, determines how much excess data to cut off due to
        # how species and time are initialized to 0
        timeStopped = 0

        if resume is not None:
            if resume[0].model != model:
                raise gillespyError.ModelError('When resuming, one must not alter the model being resumed.')
            if t < resume['time'][-1]:
                raise gillespyError.ExecutionError(
                    "'t' must be greater than previous simulations end time, or set in the run() method as the "
                    "simulations next end time")

        random.seed(seed)

        species_mappings, species, parameter_mappings, number_species = nputils.numpy_initialization(model)

        # create dictionary of all constant parameters for propensity evaluation
        parameters = {'V': model.volume}
        for paramName, param in model.listOfParameters.items():
            parameters[parameter_mappings[paramName]] = param.value

        # create mapping of reaction dictionary to array indices
        reactions = list(model.listOfReactions.keys())

        # create mapping of reactions, and which reactions depend on their reactants/products
        dependent_rxns = nputils.dependency_grapher(model, reactions)
        number_reactions = len(reactions)
        propensity_functions = {}

        # create an array mapping reactions to species modified
        species_changes = np.zeros((number_reactions, number_species))

        # pre-evaluate propensity equations from strings:
        for i, reaction in enumerate(reactions):
            # replace all references to species with array indices
            for j, spec in enumerate(species):
                species_changes[i][j] = model.listOfReactions[reaction].products.get(model.listOfSpecies[spec], 0) \
                                        - model.listOfReactions[reaction].reactants.get(model.listOfSpecies[spec], 0)
                if debug:
                    print('species_changes: {0},i={1}, j={2}... {3}'.format(species, i, j, species_changes[i][j]))
            propensity_functions[reaction] = [eval('lambda S:' + model.listOfReactions[reaction].
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

            for spec in model.listOfSpecies:
                if resume is not None:
                    curr_state[0][spec] = resume[spec][-1]
                else:
                    curr_state[0][spec] = model.listOfSpecies[spec].initial_value

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
                curr_time[0] += -math.log(random.random()) / propensity_sum
                total_time[0] += -math.log(random.random()) / propensity_sum
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

                        for i,spec in enumerate(model.listOfSpecies):
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
