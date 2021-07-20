"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from collections import OrderedDict
from gillespy2.core import Species, Reaction, Parameter, Model


class SanitizedModel:
    """
    Utility class for preparing a model to be passed to a C++ solver.
    Ensures that any sanitized mappings are consistent between different mappings.
    Wraps around an existing GillesPy2 model to provide sanitized mappings.

    :param model: GillesPy2 model to produce sanitized mappings of.
    :type model: gillespy2.Model
    """

    def __init__(self, model: Model):
        self.species: "OrderedDict[str, Species]" = OrderedDict()
        self.species_names = model.sanitized_species_names()
        for species_name, sanitized_name in self.species_names.items():
            self.species[sanitized_name] = model.get_species(species_name)

        self.parameters: "OrderedDict[str, Parameter]" = OrderedDict()
        self.parameter_names = model.sanitized_parameter_names()
        for parameter_name, sanitized_name in self.parameter_names.items():
            self.parameters[sanitized_name] = model.get_parameter(parameter_name) \
                if parameter_name != "vol" else Parameter(name="V", expression=model.volume)

        self.propensities: "OrderedDict[str, str]" = OrderedDict()
        self.ode_propensities: "OrderedDict[str, str]" = OrderedDict()
        self.reactions: "OrderedDict[str, dict[str, int]]" = OrderedDict()
        for reaction_name, reaction in model.get_all_reactions().items():
            self.propensities[reaction_name] = reaction.sanitized_propensity_function(
                self.species_names, self.parameter_names, ode=False)
            self.ode_propensities[reaction_name] = reaction.sanitized_propensity_function(
                self.species_names, self.parameter_names, ode=True)

            self.reactions[reaction_name] = {spec: int(0) for spec in self.species_names.values()}
            for reactant, stoich_value in reaction.reactants.items():
                if isinstance(reactant, Species):
                    reactant = self.species_names[reactant.name]
                self.reactions[reaction_name][reactant] -= int(stoich_value)

            for product, stoich_value in reaction.products.items():
                if isinstance(product, Species):
                    product = self.species_names[product.name]
                self.reactions[reaction_name][product] += int(stoich_value)


def write_template(path: str, model: Model, variable=False):
    """
    Write template definitions to a specified path. If the file does not exist, one is created.
    Definitions are written in `#define KEY VALUES` format.

    :param model: Model to get the macro definitions of.
    :type model: gillespy2.Model

    :param variable: Set to true to allow for non-constant parameter values.
    :type variable: bool
    """

    # Get a dictionary of model defines and transform into a list of strings in
    # `#define KEY VALUE` format.
    defines = get_model_defines(model, variable)
    write_definitions(path, defines)

def write_definitions(path: str, defines: "dict[str, str]"):
    """
    """
    # Definition dict is transformed into a list of C++ macro definitions, with:
    # `#define KEY VALUE` format.
    template_lines = [(f"#define {key} {value}\n") for key, value in defines.items()]

    # Write generated lines to the template file.
    with open(path, "w") as template_file:
        template_file.writelines(template_lines)

def get_model_defines(model: Model, variable=False) -> "dict[str, str]":
    """
    Creates a dictionary of C++ macro definitions from the given model.
    The keys of the dictionary contain the name of the macro definition.
    The values of the dictionary are the values their corresponding macro should be defined to.

    :param model: Model to get the macro definitions of.
    :type model: gillespy2.Model

    :param variable: Set to true to allow for non-constant parameter values.
    :type variable: bool

    :returns: Dictionary of fully-formatted macro definitions.
    """
    results = {}
    sanitized_model = SanitizedModel(model)

    # Get definitions for variables
    parameter_definitions = template_def_variables(sanitized_model, variable)
    results.update(parameter_definitions)

    # Get definitions for species
    species_definitions = template_def_species(sanitized_model)
    results.update(species_definitions)

    # Get definitions for reactions
    reaction_definitions = template_def_reactions(sanitized_model)
    results.update(reaction_definitions)

    # Get definitions for propensities
    stoch_propensity_definitions = template_def_propensities(sanitized_model, ode=False)
    ode_propensity_definitions = template_def_propensities(sanitized_model, ode=True)
    results.update(stoch_propensity_definitions)
    results.update(ode_propensity_definitions)

    return results


def template_def_variables(model: SanitizedModel, variable=False) -> "dict[str, str]":
    """
    Formats the relevant parameters to be passed to a C++ simulation template.
    Passed dictionaries/lists are assumed to be sanitized and sorted.

    :param model: Sanitized model containing runtime definitions.
    :type model: SanitizedModel

    :param variable: Indicates whether the model's variables should be modifiable at runtime.
    :type variable: bool

    :returns: The result as a list of tuples containing key-value pairs to be templated.
    """
    # Entries get defined as constant if variable is set to true.
    parameter_type = "VARIABLE" if variable else "CONSTANT"
    # Parameter entries, parsed and formatted
    parameter_set = []
    for param_name, parameter in model.parameters.items():
        parameter_set.append(f"{parameter_type}({param_name},{parameter.expression})")

    return {
        "GPY_PARAMETER_VALUES": " ".join(parameter_set)
    }


def template_def_species(model: SanitizedModel) -> "dict[str, str]":
    """
    Passed dictionaries/lists are assumed to be sanitized and sorted.
    Formats the relevant species data to be passed to a C++ simulation template.

    :param model: Sanitized model containing runtime definitions.
    :type model: SanitizedModel

    :returns: Dictionary of macro definitions for species and data related to species.
    """
    # Parse and format species initial populations
    num_species = len(model.species)
    populations = OrderedDict()

    for spec_name, spec in model.species.items():
        populations[spec_name] = str(spec.initial_value)
    # Species names, parsed and formatted
    sanitized_names = [f"SPECIES_NAME({name})" for name in populations.keys()]
    populations = f"{{{','.join(populations.values())}}}"

    # Match each parameter with its macro definition name
    return {
        "GPY_INIT_POPULATIONS": populations,
        "GPY_NUM_SPECIES": num_species,
        "GPY_SPECIES_NAMES": " ".join(sanitized_names)
    }


def template_def_reactions(model: SanitizedModel, ode=False) -> "dict[str, str]":
    """
    Passed dictionaries/lists are assumed to be sanitized and sorted.
    Formats the relevant reactions and propensities to be passed to a C++ simulation template.

    :param model: Sanitized model containing runtime definitions.
    :type model: SanitizedModel

    :returns: Dictionary of macro definitions for reactions.
    """
    num_reactions = str(len(model.reactions))
    reaction_set = OrderedDict()

    for rxn_name, reaction in model.reactions.items():
        stoich = [str(int(reaction[species])) for species in model.species_names.values()]
        reaction_set[rxn_name] = f"{{{','.join(stoich)}}}"

    reaction_names = " ".join([f"REACTION_NAME({rxn})" for rxn in reaction_set.keys()])
    reaction_set = f"{{{','.join(reaction_set.values())}}}"

    return {
        "GPY_NUM_REACTIONS": num_reactions,
        "GPY_REACTIONS": reaction_set,
        "GPY_REACTION_NAMES": reaction_names,
    }


def template_def_propensities(model: SanitizedModel, ode=False) -> "dict[str, str]":
    """
    Formats the given list of pre-sorted, pre-sanitized propensities.

    :param model: Sanitized model containing runtime definitions.
    :type model: SanitizedModel

    :returns: Dictionary containing propensity macro definitions.
    """
    def_keyword = "ODE_PROPENSITIES" if ode else "PROPENSITIES"
    propensities = model.ode_propensities if ode else model.propensities
    propensities = [f"PROPENSITY({rxn_i},{propensities[rxn]})" for rxn_i, rxn in enumerate(model.reactions.keys())]

    return {
        f"GPY_{def_keyword}": " ".join(propensities)
    }
