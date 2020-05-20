
import os #for getting directories for C++ files
import shutil #for deleting/copying files
import numpy as np
from gillespy2.core.gillespyError import ExecutionError



"""
This file contains various functions used in the ssa_c_solver, variable_ssa_c_solver, and numpy solvers.


C SOLVER FUNCTIONS BELOW
"""


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


def _parse_binary_output(results_buffer, number_of_trajectories, number_timesteps, number_species):
    trajectory_base = np.empty((number_of_trajectories, number_timesteps, number_species+1))
    step_size = number_species * number_of_trajectories + 1 #1 for timestep
    data = np.frombuffer(results_buffer, dtype=np.float64)
    assert(len(data) == (number_of_trajectories*number_timesteps*number_species + number_timesteps))
    for timestep in range(number_timesteps):
        index = step_size * timestep
        trajectory_base[:, timestep, 0] = data[index]
        index += 1
        for trajectory in range(number_of_trajectories):
            for species in range(number_species):
                trajectory_base[trajectory, timestep, 1 + species] = data[index + species]
            index += number_species
    return trajectory_base

def c_solver_results(return_code,stdout,number_of_trajectories,number_timesteps,model,show_labels):
    if return_code in [0, 33]:
        trajectory_base = _parse_binary_output(stdout, number_of_trajectories, number_timesteps, len(model.species))
        # Format results
        if show_labels:
            model.simulation_data = []
            for trajectory in range(number_of_trajectories):
                data = {'time': trajectory_base[trajectory, :, 0]}
                for i in range(len(model.species)):
                    data[model.species[i]] = trajectory_base[trajectory, :, i + 1]
                model.simulation_data.append(data)
        else:
            model.simulation_data = trajectory_base
    else:
        raise ExecutionError(
            "Error encountered while running simulation C++ file:\nReturn code: {0}.\nError:\n{1}\n".format(
                model.simulation.returncode, model.simulation.stderr))
    return model.simulation_data, return_code


"""
NUMPY SOLVER FUNCTIONS BELOW
"""

def numpyresults(data, species, number_species, trajectory,simulation_data):
    for i in range(number_species):
        data[species[i]] = trajectory[:, i + 1]
    simulation_data.append(data)
    return simulation_data