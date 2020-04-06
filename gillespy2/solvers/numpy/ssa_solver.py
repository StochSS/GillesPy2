from threading import Thread, Event
from gillespy2.core import GillesPySolver, Model, Reaction, log
import random
import math
import numpy as np


class NumPySSASolver(GillesPySolver):
    name = "NumPySSASolver"
    rc = 0
    stop_event = None
    result = None

    def __init__(self):
        name = 'NumPySSASolver'
        rc = 0
        stop_event = None
        result = None

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05,
                        seed=None, debug=False, show_labels=True, timeout=None, **kwargs):
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
        :param show_labels: Use names of species as index of result object rather than position numbers.
        :return: a list of each trajectory simulated.
        """


        if not isinstance(self, NumPySSASolver):
            self = NumPySSASolver()

        self.stop_event = Event()
        if timeout is not None and timeout <= 0: timeout = None

        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))
        sim_thread = Thread(target=self.__run, args=(model,), kwargs={'t':t,
                                        'number_of_trajectories':number_of_trajectories,
                                        'increment':increment, 'seed':seed,
                                        'debug':debug, 'show_labels':show_labels,
                                        'timeout':timeout})
        sim_thread.start()
        sim_thread.join(timeout=timeout)
        self.stop_event.set()
        while self.result is None: pass
        return self.result, self.rc

    def __run(self, model, t=20, number_of_trajectories=1, increment=0.05,
                    seed=None, debug=False, show_labels=True, timeout=None):

        random.seed(seed)
        # create mapping of species dictionary to array indices
        species_mappings = model.sanitized_species_names()
        species = list(species_mappings.keys())
        parameter_mappings = model.sanitized_parameter_names()
        number_species = len(species)

        # create numpy array for timeline
        timeline = np.linspace(0, t, int(round(t / increment + 1)))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline
        # copy initial populations to base

        for i, s in enumerate(species):
            trajectory_base[:, 0, i + 1] = model.listOfSpecies[s].initial_value
            # create dictionary of all constant parameters for propensity evaluation
        parameters = {'V': model.volume}
        for paramName, param in model.listOfParameters.items():
            parameters[parameter_mappings[paramName]] = param.value

        # create mapping of reaction dictionary to array indices
        reactions = list(model.listOfReactions.keys())
        number_reactions = len(reactions)
        propensity_functions = []
        # create an array mapping reactions to species modified
        species_changes = np.zeros((number_reactions, number_species))
        # pre-evaluate propensity equations from strings:
        for i, reaction in enumerate(reactions):
            # replace all references to species with array indices
            for j, spec in enumerate(species):
                species_changes[i][j] = model.listOfReactions[reaction].products.get(model.listOfSpecies[spec], 0) - model.listOfReactions[reaction].reactants.get(model.listOfSpecies[spec], 0)
                if debug:
                    print('species_changes: {0},i={1}, j={2}... {3}'.format(species, i, j, species_changes[i][j]))
            propensity_functions.append(eval('lambda S:' + model.listOfReactions[reaction].sanitized_propensity_function(species_mappings, parameter_mappings), parameters))
        if debug:
            print('propensity_functions', propensity_functions)
        # begin simulating each trajectory
        simulation_data = []
        for trajectory_num in range(number_of_trajectories):
            if self.stop_event.is_set():
                self.rc = 33
                break
            # copy initial state data
            trajectory = trajectory_base[trajectory_num]
            entry_count = 1
            current_time = 0
            current_state = np.copy(trajectory[0, 1:])
            propensity_sums = np.zeros(number_reactions)
            # calculate initial propensity sums
            while entry_count < timeline.size:
                if self.stop_event.is_set(): 
                    self.rc = 33
                    break
                # determine next reaction
                for i in range(number_reactions):
                    propensity_sums[i] = propensity_functions[i](current_state)
                    if debug:
                        print('propensity: ', propensity_sums[i])
                propensity_sum = np.sum(propensity_sums)
                if debug:
                    print('propensity_sum: ', propensity_sum)
                # if no more reactions, quit
                if propensity_sum <= 0:
                    trajectory[entry_count:, 1:] = current_state
                    break
                cumulative_sum = random.uniform(0, propensity_sum)
                current_time += -math.log(random.random()) / propensity_sum
                if debug:
                    print('cumulative sum: ', cumulative_sum)
                    print('entry count: ', entry_count)
                    print('timeline.size: ', timeline.size)
                    print('current_time: ', current_time)
                # determine time passed in this reaction
                while entry_count < timeline.size and timeline[entry_count] <= current_time:
                    if self.stop_event.is_set(): 
                        self.rc = 33
                        break
                    trajectory[entry_count, 1:] = current_state
                    entry_count += 1
                for potential_reaction in range(number_reactions):
                    cumulative_sum -= propensity_sums[potential_reaction]
                    if debug:
                        print('if <=0, fire: ', cumulative_sum)
                    if cumulative_sum <= 0:
                        current_state += species_changes[potential_reaction]
                        if debug:
                            print('current state: ', current_state)
                            print('species_changes: ', species_changes)
                            print('updating: ', potential_reaction)
                        # recompute propensities as needed
                        for i in range(number_reactions):
                            propensity_sums[i] = propensity_functions[i](current_state)
                            if debug:
                                print('new propensity sum: ', propensity_sums[i])
                        break
            if show_labels:
                data = {
                    'time': timeline
                }
                for i in range(number_species):
                    data[species[i]] = trajectory[:, i+1]
                simulation_data.append(data)
            else:
                simulation_data = trajectory_base
        self.result = simulation_data
        return self.result, self.rc
