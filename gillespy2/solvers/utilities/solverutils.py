import os  # for getting directories for C++ files
import shutil  # for deleting/copying files
import ast  # for dependency graphing
import numpy as np
from gillespy2.core import log, Species


"""
This file contains various functions used in the ssa_c_solver, variable_ssa_c_solver, and numpy solvers.


C SOLVER FUNCTIONS BELOW
"""


def find_time(array, value):
    """
    Finds the index of the closest value in the array parameter, to the value parameter
    :param array: results['time'] array, input to find index of closest 'value' parameter
    :type array: numpy.ndarray
    :param value: Value in which to find the closest values index in the array parameter.
    :type value: float
    :return: Integer index, the index of the closest value to 'value' parameter.
    """
    index = np.searchsorted(array, value, side="left")
    return index


def _copy_files(destination, GILLESPY_C_DIRECTORY):
    src_files = os.listdir(GILLESPY_C_DIRECTORY)
    for src_file in src_files:
        src_file = os.path.join(GILLESPY_C_DIRECTORY, src_file)
        if os.path.isfile(src_file):
            shutil.copy(src_file, destination)


def _write_propensity(outfile, model, species_mappings, parameter_mappings, reactions):
    """
    This functions writes a models propensity functions to a cpp user simulation template, for the SSACSolvers.
    :param outfile: File where the propensity function will be written to
    :param model: Model used to access species, reactions
    :param species_mappings: Sanitized species names
    :param parameter_mappings: Sanitized parameter names
    :param reactions: Names of reactions
    :type reactions: str
    """
    for i in range(len(reactions)):
        # Write switch statement case for reaction
        outfile.write("""
        case {0}:
            return {1};
        """.format(i, model.listOfReactions[reactions[i]].sanitized_propensity_function(species_mappings,
                                                                                        parameter_mappings)))


def _write_reactions(outfile, model, reactions, species):
    customrxns = {}
    for i in range(len(reactions)):
        reaction = model.listOfReactions[reactions[i]]
        if reaction.type == 'customized':
            customrxns[i] = species_parse(model, reaction.propensity_function)
        for j in range(len(species)):
            change = (reaction.products.get(model.listOfSpecies[species[j]], 0)) - (reaction.reactants.get
                                                                                    (model.listOfSpecies[species[j]], 0)
                                                                                    )
            if change != 0:
                outfile.write("model.reactions[{0}].species_change[{1}] = {2};\n".format(i, j, change))

    for i in customrxns.keys():
        for j in range(len(reactions)):
            if i == j:
                continue
            if any(elem in customrxns[i] for elem in list(model.listOfReactions[reactions[j]].reactants)) or \
                    any(elem in customrxns[i] for elem in list(model.listOfReactions[reactions[j]].products)):
                outfile.write("model.reactions[{0}].affected_reactions.push_back({1});\n".format(i, j))


def _parse_binary_output(results_buffer, number_of_trajectories, number_timesteps, number_species, pause=False):
    """
    This function reads binary output from a CPP simulation
    :param results_buffer: stdout of the CPP simulation ran
    :param number_of_trajectories: Total number of trajectories for a simulation
    :type number_of_trajectories: int
    :param number_timesteps: How many steps for a given simulation
    :type number_timesteps: int
    :param number_species: Total number of species in a model
    :param pause: Whether or not a model was paused, set to true when simulation was sent a KeyBoardInterrupt or
    timeout.
    :return: Trajectories for a simulation, and time that simulation was stopped, if sent a keyboardinterrupt or
    timeout.
    """
    trajectory_base = np.empty((number_of_trajectories, number_timesteps, number_species+1))
    step_size = number_species * number_of_trajectories + 1  #1 for timestep
    data = np.frombuffer(results_buffer, dtype=np.float64)
    # Timestopped is added to the end of the data, when a simulation completes or is paused
    if pause:
        timeStopped = data[-1]
    else:
        timeStopped = 0
    assert(len(data) == (number_of_trajectories*number_timesteps*number_species + number_timesteps)+1)
    for timestep in range(number_timesteps):
        index = step_size * timestep
        trajectory_base[:, timestep, 0] = data[index]
        index += 1
        for trajectory in range(number_of_trajectories):
            for species in range(number_species):
                trajectory_base[trajectory, timestep, 1 + species] = data[index + species]
            index += number_species
    return trajectory_base, timeStopped


def c_solver_resume(timeStopped, simulation_data, t, resume=None):
    """
    If a simulation is being resumed from a previous simulation, this function is called in the VariableSSACSolver,
    or SSACSolver
    :param timeStopped: The time that a simulation was stopped, originally attained from the results_buffer returned
    by the CPP simulation.
    :param simulation_data: The current simulation data, attained after parsing the results in the VariableSSACSolver or
    SSACSolver.
    :param t: The end time for the resume simulation, originally set in model.run(t=...)
    :param resume: The previous simulations data
    :type resume: gillespy2.core.result object
    :return: Combined data of the previous simulation, and the current simulation
    """

    # If simulation was paused/KeyboardInterrupt
    if timeStopped != 0:
        cutoff = find_time(simulation_data[0]['time'], timeStopped)
        if cutoff == 0 or cutoff == 1:
            log.warning('You have paused the simulation too early, and no points have been calculated past'
                        ' initial values. A graphic display will not produce expected results.')
        else:
            cutoff -= 1
        for i in simulation_data[0]:
            simulation_data[0][i] = simulation_data[0][i][:cutoff]

    if resume is not None:
        resumeTime = float(resume['time'][-1])
        step = resumeTime - resume['time'][-2]
        if timeStopped == 0:
            timeSpan = np.arange(resumeTime, t + resumeTime + step, step)
        else:
            timeSpan = np.arange(resumeTime + step, timeStopped + resumeTime + step, step)
        simulation_data[0]['time'] = timeSpan

    if resume is not None:
        # If resuming, combine old pause with new data, and delete any excess null data
        for i in simulation_data[0]:
            oldData = resume[i]
            newData = simulation_data[0][i]
            simulation_data[0][i] = np.concatenate((oldData, newData), axis=None)
            if len(simulation_data[0]['time']) != len(simulation_data[0][i]):
                simulation_data[0]['time'] = simulation_data[0]['time'][:-1]
    return simulation_data


"""
NUMPY SOLVER UTILITIES BELOW
"""


def numpy_initialization(model):
    species_mappings = model.sanitized_species_names()
    species = list(species_mappings.keys())
    parameter_mappings = model.sanitized_parameter_names()
    number_species = len(species)
    return species_mappings, species, parameter_mappings, number_species


def numpy_trajectory_base_initialization(model, number_of_trajectories, timeline, species, resume=None):
    trajectory_base = np.zeros((number_of_trajectories, timeline.size, len(species) + 1))

    # copy time values to all trajectory row starts
    trajectory_base[:, :, 0] = timeline
    tmpSpecies = {}
    # copy initial populations to base
    if resume is not None:
        # Set initial values of species to where last left off
        for i in species:
            tmpSpecies[i] = resume[i][-1]
        for i, s in enumerate(species):
            trajectory_base[:, 0, i + 1] = tmpSpecies[s]
    else:
        for i, s in enumerate(species):
            trajectory_base[:, 0, i + 1] = model.listOfSpecies[s].initial_value

    return trajectory_base, tmpSpecies


def numpy_resume(timeStopped, simulation_data, resume=None):
    """
    Helper function for when resuming a simulation in a numpy based solver
    :param timeStopped: The time in which the simulation was stopped.
    :param simulation_data: The current models simulation data, after being parsed in the numpy solver of choice.
    :param resume: The previous simulations data, that is being resumed
    :type resume: gillespy2.core.results object
    :return: Combined simulation data, the old resume data and the current simulation data.
    """
    if timeStopped != 0:
        if timeStopped != simulation_data[0]['time'][-1]:
            tester = np.where(simulation_data[0]['time'] > timeStopped)[0].size
            index = np.where(simulation_data[0]['time'] == timeStopped)[0][0]
        if tester > 0:
            for i in simulation_data[0]:
                simulation_data[0][i] = simulation_data[0][i][:index]
    if resume is not None:
        # If resuming, combine old pause with new data, and delete any excess null data
        for i in simulation_data[0]:
            oldData = resume[i][:-1]
            newData = simulation_data[0][i]
            simulation_data[0][i] = np.concatenate((oldData, newData), axis=None)

    return simulation_data


"""
Below are two functions used for creating dependency graphs in the C solvers, and Numpy Solvers.
"""


def species_parse(model, custom_prop_fun):
    """
    This function uses Pythons AST module to parse custom propensity function, looking for Species in a model
    :param model: Model to be checked for species
    :param custom_prop_fun: The custom propensity function to be parsed
    :return: List of species objects that are found in a custom propensity function
    """
    parsed_species = []

    class SpeciesParser(ast.NodeTransformer):
        def visit_Name(self, node):
            if isinstance(model.get_element(node.id), Species):
                parsed_species.append(model.get_element(node.id))

    expr = custom_prop_fun
    expr = ast.parse(expr, mode='eval')
    expr = SpeciesParser().visit(expr)
    return parsed_species


def dependency_grapher(model, reactions):
    """
    This function returns a dependency graph for a models reactions in the form of a
    dictionary containing {species name: {'dependencies'}:[list of reaction names]}.
    :param model: Model to used to create a reaction dependency graph
    :param reactions: list(model.listOfReactions)
    :return: Dependency graph dictionary
    """
    dependent_rxns = {}
    for i in reactions:
        cust_spec = []
        if model.listOfReactions[i].type == 'customized':
            cust_spec = (species_parse(model, model.listOfReactions[i].propensity_function))

        for j in reactions:

            if i not in dependent_rxns:
                dependent_rxns[i] = {'dependencies': []}
            if j not in dependent_rxns:
                dependent_rxns[j] = {'dependencies': []}
            if i == j:
                continue

            reactantsI = list(model.listOfReactions[i].reactants.keys())
            reactantsJ = list(model.listOfReactions[j].reactants.keys())

            if j not in dependent_rxns[i]['dependencies']:
                if any(elem in reactantsI for elem in reactantsJ):
                    if i not in dependent_rxns[j]['dependencies']:
                        dependent_rxns[j]['dependencies'].append(i)
                    dependent_rxns[i]['dependencies'].append(j)

            if i not in dependent_rxns[j]['dependencies']:
                if any(elem in list(model.listOfReactions[i].products.keys()) for elem in
                       list(model.listOfReactions[j].reactants.keys())):
                    dependent_rxns[j]['dependencies'].append(i)

            if cust_spec:
                if any(elem in cust_spec for elem in list(model.listOfReactions[j].reactants)) or any \
                            (elem in cust_spec for elem in list(model.listOfReactions[j].products)):
                    dependent_rxns[i]['dependencies'].append(j)

    return dependent_rxns

