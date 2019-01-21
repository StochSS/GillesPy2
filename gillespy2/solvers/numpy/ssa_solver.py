from gillespy2.core import GillesPySolver, Model, Reaction
import random
import math
import numpy as np


class NumPySSASolver(GillesPySolver):
    name = "NumPySSASolver"

    @staticmethod
    def run(model, t=20, number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=False, **kwargs):
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
        # create mapping of species dictionary to array indices
        species_mappings = model.sanitized_species_names()
        species = list(species_mappings.keys())
        parameter_mappings = model.sanitized_parameter_names()
        number_species = len(species)

        # create numpy array for timeline
        timeline = np.linspace(0, t, (t // increment + 1))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.empty((number_of_trajectories, timeline.size, number_species + 1))

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
            propensity_functions.append(eval('lambda S:' + model.listOfReactions[reaction].sanitized_propensity_function(species_mappings, parameter_mappings), parameters))
        # begin simulating each trajectory
        simulation_data = []
        for trajectory_num in range(number_of_trajectories):
            # copy initial state data
            trajectory = trajectory_base[trajectory_num]
            entry_count = 1
            current_time = 0
            current_state = np.copy(trajectory[0, 1:])
            propensity_sums = np.zeros(number_reactions)
            # calculate initial propensity sums
            while entry_count < timeline.size:
                # determine next reaction
                for i in range(number_reactions):
                    propensity_sums[i] = propensity_functions[i](current_state)
                propensity_sum = np.sum(propensity_sums)
                # if no more reactions, quit
                if propensity_sum <= 0:
                    trajectory[entry_count:, 1:] = current_state
                    break
                cumulative_sum = random.uniform(0, propensity_sum)
                current_time += -math.log(random.random()) / propensity_sum
                # determine time passed in this reaction
                while entry_count < timeline.size and timeline[entry_count] <= current_time:
                    trajectory[entry_count, 1:] = current_state
                    entry_count += 1
                for potential_reaction in range(number_reactions):
                    cumulative_sum -= propensity_sums[potential_reaction]
                    if cumulative_sum <= 0:
                        current_state += species_changes[potential_reaction]
                        # recompute propensities as needed
                        for i in range(number_reactions):
                            propensity_sums[i] = propensity_functions[i](current_state)
                        break
            if show_labels:
                data = {
                    'time': timeline
                }
                for i in range(number_species):
                    data[species[i]] = trajectory[:, i+1]
                simulation_data.append(data)
            else:
                simulation_data.append(trajectory)
        return simulation_data
