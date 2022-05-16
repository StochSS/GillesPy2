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

import ast  # for dependency graphing
import numpy as np
from gillespy2.core import log, Species
from gillespy2.core import ModelError

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
    Helper function for when resuming a simulation in a numpy based solver.

    :param timeStopped: The time in which the simulation was stopped.
    :param simulation_data: The current models simulation data, after being parsed in the numpy solver of choice.

    :param resume: The previous simulations data, that is being resumed
    :type resume: gillespy2.core.Results

    :returns: Combined simulation data, the old resume data and the current simulation data.
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
VARIABLE SOLVER METHODS
"""


def update_species_init_values(listOfSpecies, species, variables, resume = None):
    # Update Species Initial Values
    populations = ''
    for i in range(len(species) - 1):
        if species[i] in variables:
            populations += '{} '.format(float(variables[species[i]]))
        else:
            if resume is not None:
                populations += '{} '.format(float(resume[species[i]][-1]))
            else:
                populations += '{} '.format(float(listOfSpecies[species[i]].initial_value))
    if species[-1] in variables:
        populations += '{}'.format(float(variables[species[-1]]))
    else:
        if resume is not None:
            populations += '{} '.format(float(resume[species[-1]][-1]))
        else:
            populations += '{}'.format(float(listOfSpecies[species[-1]].initial_value))
    return populations

def change_param_values(listOfParameters, parameters, volume, variables):
    # Update Parameter Values
    parameter_values = ''
    for i in range(len(parameters) - 1):
        if parameters[i] in variables:
            parameter_values += '{} '.format(variables[parameters[i]])
        else:
            if parameters[i] == 'vol':
                parameter_values += '{} '.format(volume)
            else:
                parameter_values += '{} '.format(listOfParameters[parameters[i]].expression)
    if parameters[-1] in variables:
        parameter_values += '{}'.format(variables[parameters[-1]])
    else:
        if parameters[-1] == 'vol':
            parameter_values += '{}'.format(volume)
        else:
            parameter_values += '{}'.format(listOfParameters[parameters[-1]].expression)
    return parameter_values

"""
Below are two functions used for creating dependency graphs in the C solvers, and Numpy Solvers.
"""

def species_parse(model, custom_prop_fun):
    """
    This function uses Pythons AST module to parse custom propensity function, looking for Species in a model

    :param model: Model to be checked for species
    :param custom_prop_fun: The custom propensity function to be parsed

    :returns: List of species objects that are found in a custom propensity function
    """
    parsed_species = []

    class SpeciesParser(ast.NodeTransformer):
        def visit_Name(self, node):
            try:
                if isinstance(model.get_element(node.id), Species):
                    parsed_species.append(model.get_element(node.id))
            except ModelError:
                pass

    expr = custom_prop_fun
    expr = ast.parse(expr, mode='eval')
    expr = SpeciesParser().visit(expr)
    return parsed_species


def dependency_grapher(model, reactions):
    """
    This function returns a dependency graph for a models reactions in the form of a
    dictionary containing {species name: {'dependencies'}:[list of reaction names]}.

    :param model: Model to used to create a reaction dependency graph
    :param reactions: list[model.listOfReactions]

    :returns: Dependency graph dictionary
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

