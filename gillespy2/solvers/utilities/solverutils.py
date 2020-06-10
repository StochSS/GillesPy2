import os #for getting directories for C++ files
import shutil #for deleting/copying files
import numpy as np
from gillespy2.core import log



"""
This file contains various functions used in the ssa_c_solver, variable_ssa_c_solver, and numpy solvers.


C SOLVER FUNCTIONS BELOW
"""
def find_time(array,value):
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

def _copy_files(destination,GILLESPY_C_DIRECTORY):
    src_files = os.listdir(GILLESPY_C_DIRECTORY)
    for src_file in src_files:
        src_file = os.path.join(GILLESPY_C_DIRECTORY, src_file)
        if os.path.isfile(src_file):
            shutil.copy(src_file, destination)

def _write_propensity(outfile, model, species_mappings, parameter_mappings, reactions):
    for i in range(len(reactions)):
        # Write switch statement case for reaction
        outfile.write("""
        case {0}:
            return {1};
        """.format(i, model.listOfReactions[reactions[i]].sanitized_propensity_function(species_mappings, parameter_mappings)))


def _write_reactions(outfile, model, reactions, species):
    for i in range(len(reactions)):
        reaction = model.listOfReactions[reactions[i]]
        for j in range(len(species)):
            change = (reaction.products.get(model.listOfSpecies[species[j]], 0)) - (reaction.reactants.get(model.listOfSpecies[species[j]], 0))
            if change != 0:
                outfile.write("model.reactions[{0}].species_change[{1}] = {2};\n".format(i, j, change))


def _parse_output(results, number_of_trajectories, number_timesteps, number_species):
    trajectory_base = np.empty((number_of_trajectories, number_timesteps, number_species+1))
    for timestep in range(number_timesteps):
        values = results[timestep].split(" ")
        trajectory_base[:, timestep, 0] = float(values[0])
        index = 1
        for trajectory in range(number_of_trajectories):
            for species in range(number_species):
                trajectory_base[trajectory, timestep, 1 + species] = float(values[index+species])
            index += number_species
    return trajectory_base


def _parse_binary_output(results_buffer, number_of_trajectories, number_timesteps, number_species,pause=False):
    trajectory_base = np.empty((number_of_trajectories, number_timesteps, number_species+1))
    step_size = number_species * number_of_trajectories + 1 #1 for timestep
    data = np.frombuffer(results_buffer, dtype=np.float64)
    #Timestopped is added to the end of the data, when a simulation completes or is paused
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
