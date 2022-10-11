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

from collections import OrderedDict
from typing import Optional
from gillespy2.core import Species, Reaction, Parameter, Model, RateRule
from gillespy2.solvers.cpp.build.expression import Expression
from gillespy2.core import log
from gillespy2.core.gillespyError import SimulationError

import math


class SanitizedModel:
    """
    Utility class for preparing a model to be passed to a C++ solver.
    Ensures that any sanitized mappings are consistent between different mappings.
    Wraps around an existing GillesPy2 model to provide sanitized mappings.

    :param model: GillesPy2 model to produce sanitized mappings of.
    :type model: gillespy2.Model
    """
    reserved_names = {
        "t": "t",
    }

    # Global functions that aren't present in the `math` package,
    # as well as functions in Python that have a different name in C++.
    function_map = {
        "abs": "abs",
        "round": "round",
    }

    def __init__(self, model: Model, variable=False):
        self.model = model
        self.variable = variable

        self.species: "OrderedDict[str, Species]" = OrderedDict()
        self.species_names = model.sanitized_species_names()
        self.species_id: "OrderedDict[str, int]" = OrderedDict()
        for spec_id, spec_entry in enumerate(self.species_names.items()):
            species_name, sanitized_name = spec_entry
            self.species[sanitized_name] = model.get_species(species_name)
            self.species_id[species_name] = spec_id

        self.parameters: "OrderedDict[str, Parameter]" = OrderedDict()
        self.parameter_names: "OrderedDict[str, str]" = OrderedDict()
        self.parameter_names["vol"] = "P[0]" if variable else "C[0]"
        self.parameter_id: "OrderedDict[str, int]" = OrderedDict()
        for param_id, param_name in enumerate(model.listOfParameters.keys(), start=1):
            if param_name not in self.parameter_names:
                self.parameter_names[param_name] = f"P[{param_id}]" if variable else f"C[{param_id}]"
                self.parameter_id[param_name] = param_id
        for parameter_name, sanitized_name in self.parameter_names.items():
            self.parameters[sanitized_name] = model.get_parameter(parameter_name) \
                if parameter_name != "vol" else Parameter(name=sanitized_name, expression=str(model.volume))

        base_namespace = {
            # ORDER IS IMPORTANT HERE!
            # All "system" namespace entries should always be first.
            # Otherwise, user-defined identifiers (like, for example, "gamma") might get overwritten.
            **{name: name for name in math.__dict__.keys()},
            **self.function_map,
            **self.species_names,
            **self.parameter_names,
            **self.reserved_names,
        }
        self.expr = Expression(namespace=base_namespace, blacklist=["="], sanitize=True)

        # SSA Propensities: Maps reaction names to their corresponding propensity function.
        self.propensities: "OrderedDict[str, str]" = OrderedDict()
        # ODE Propensities: Maps reaction names to their corresponding mass-action rate expression.
        self.ode_propensities: "OrderedDict[str, str]" = OrderedDict()
        # Reactions: maps reaction names to their stoichiometry matrix.
        # Stoichiometry matrix maps a sanitized species name to its stoichiometry.
        self.reactions: "OrderedDict[str, dict[str, int]]" = OrderedDict()
        self.reaction_reactants: "OrderedDict[str, dict[str, int]]" = OrderedDict()
        self.reaction_products: "OrderedDict[str, dict[str, int]]" = OrderedDict()
        # Rate Rules: maps sanitized species names to their corresponding rate rule expression.
        self.rate_rules: "OrderedDict[str, str]" = OrderedDict()
        # Options: custom definitions that can be supplied by the solver, maps macros to their definitions.
        # The solver itself may use `options` to supply their own solver-specific definitions.
        self.options: "OrderedDict[str, str]" = OrderedDict()

        for reaction in model.get_all_reactions().values():
            self.use_reaction(reaction)
            self.use_propensity(reaction, ode=False)
            self.use_propensity(reaction, ode=True)

    def use_propensity(self, reaction: "Reaction", ode=False):
        """
        Populates the given reaction's propensities into the sanitized model.
        Expression conversion and sanitization are automatically applied.

        :param reaction: Reaction to generate propensities from.
        :type reaction: gillespy2.Reaction

        :param ode: Determines whether the stochastic propensity or ODE mass-action rate is used.
        If set to true, the mass-action rate is used instead of the propensity.
        """
        propensities = self.ode_propensities if ode else self.propensities
        propensity = reaction.ode_propensity_function if ode else reaction.propensity_function
        propensities[reaction.name] = self.expr.getexpr_cpp(propensity)

    def use_reaction(self, reaction: "Reaction") -> "SanitizedModel":
        """
        Adds the given reaction to the sanitized model.
        Populates the name and stoichiometry matrix into the model.

        :param reaction: Reaction to add to the sanitized model.
        :type reaction: gillespy2.Reaction
        """
        self.reactions[reaction.name] = {spec: int(0) for spec in self.species_names.values()}
        self.reaction_reactants[reaction.name] = {spec: int(0) for spec in self.species_names.values()}
        self.reaction_products[reaction.name] = {spec: int(0) for spec in self.species_names.values()}
        for reactant, stoich_value in reaction.reactants.items():
            reactant = self.species_names[reactant.name]
            self.reactions[reaction.name][reactant] -= int(stoich_value)
            self.reaction_reactants[reaction.name][reactant] = int(stoich_value)

        for product, stoich_value in reaction.products.items():
            product = self.species_names[product.name]
            self.reactions[reaction.name][product] += int(stoich_value)
            self.reaction_products[reaction.name][product] = int(stoich_value)

        return self

    def use_rate_rule(self, rate_rule: "RateRule") -> "SanitizedModel":
        """
        Attach the given rate rule to the sanitized model.
        The rate rule will automatically be validated and sanitized before being applied.

        :param rate_rule: GillesPy2 RateRule object to attach to the sanitized model.
        :type rate_rule: gillespy2.RateRule

        :returns: Pass-through of sanitized model object.
        :rtype: SanitizedModel
        """

        if isinstance(rate_rule.variable, Species):
            variable = rate_rule.variable
        elif not isinstance(rate_rule.variable, Parameter) and \
                    rate_rule.variable in self.model.listOfSpecies.keys():
            variable = self.model.get_species(rate_rule.variable)
        else:
            errmsg = """
            Parameters are not valid variables for the TauHybridCSolver.

            In order to use this variable it will need to be a gillespy2.Species.
            """
            raise SimulationError(errmsg)

        if variable.name in self.species_names:
            sanitized_name = self.species_names.get(variable.name)
            if sanitized_name in self.rate_rules:
                log.warning(f"Duplicate rate rule variable found in C++ solver: {variable}")
            rr_sanitized = self.expr.getexpr_cpp(rate_rule.formula)
            if rr_sanitized is not None:
                self.rate_rules[sanitized_name] = rr_sanitized
            else:
                log.warning(f"Could not sanitize rate rule formula expression: {rate_rule.formula}")
        return self

    def get_template(self) -> "dict[str, str]":
        """
        Creates a dictionary of C++ macro definitions from the given model.
        The keys of the dictionary contain the name of the macro definition.
        The values of the dictionary are the values their corresponding macro should be defined to.

        :param variable: Set to true to allow for non-constant parameter values.
        :type variable: bool

        :returns: Dictionary of fully-formatted macro definitions.
        :rtype: dict[str, str]
        """
        results = dict({})

        # Get definitions for variables
        parameter_definitions = template_def_variables(self, self.variable)
        results.update(parameter_definitions)

        # Get definitions for species
        species_definitions = template_def_species(self)
        results.update(species_definitions)

        # Get definitions for reactions
        reaction_definitions = template_def_reactions(self)
        results.update(reaction_definitions)

        # Get definitions for propensities
        stoch_propensity_definitions = template_def_propensities(self, ode=False)
        ode_propensity_definitions = template_def_propensities(self, ode=True)
        results.update(stoch_propensity_definitions)
        results.update(ode_propensity_definitions)

        return results

    def get_options(self) -> "Optional[dict[str, str]]":
        """
        Creates a dictionary of C++ macro definitions for optional parameters of the model.
        The keys of the dictionary contain the name of the macro definition.
        The values of the dictionary are the values their corresponding macro should be defined to.

        :returns: Dictionary of fully-formatted macro definitions.
        :rtype: dict[str, str]
        """
        options = update_model_options(self, self.options.copy())
        return options if len(options) > 0 else None


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
    Write the given key-value pairs as a C/C++ template file.
    Contents of the given `defines` dict are written as `#define {key} {value}` format.

    :param path: Absolute filepath of the header file to create (including the `.h` extension.
    :type path: str

    :param defines: Dictionary of key-value pairs representing macro definitions to create.
    :type defines: dict[str, str]
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


def update_model_options(model: "SanitizedModel", definitions: "dict[str, str]" = None) -> "dict[str, str]":
    """
    Generate the solver-specific options of the `SanitizedModel` with the given definitions.
    This includes both solver-defined custom definitions, as well as SBML features.

    :param model: Sanitized model to populate the options into.
    The model itself is not modified.
    :type model: SanitizedModel

    :param definitions: Dictionary of macro definition key-value pairs to pass to the model.
    :type definitions: dict[str, str]

    :returns: Updated dictionary of macro definitions.
    Includes both solver-specific definitions provided and SBML features on the model.
    :rtype: dict[str, str]
    """
    if definitions is None:
        definitions = {}

    if len(model.rate_rules) > 0:
        definitions.update(template_def_rate_rules(model))

    return definitions


def template_def_rate_rules(model: SanitizedModel) -> "dict[str, str]":
    """
    Generates template definitions for SBML rate rules.
    """
    rr_set = []
    for spec_i, species in enumerate(model.species_names.values()):
        if species in model.rate_rules:
            spec_rr = model.rate_rules.get(species)
            rr_set.append(f"RATE_RULE({spec_i}, {spec_rr})")

    return {
        "GPY_RATE_RULES": " ".join(rr_set)
    }

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
    for param_id, parameter in enumerate(model.parameters.values()):
        parameter_set.append(f"{parameter_type}({param_id},{parameter.expression})")

    return {
        "GPY_PARAMETER_VALUES": " ".join(parameter_set),
        # Currently assumes all variable or all constant.
        # For partially variable models, modify to compute these two separately.
        "GPY_PARAMETER_NUM_VARIABLES": str(len(parameter_set)) if variable else "0",
        "GPY_PARAMETER_NUM_CONSTANTS": str(len(parameter_set)) if not variable else "0",
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
        populations[spec_name] = str(float(spec.initial_value))
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
    reactants_set = OrderedDict()
    products_set = OrderedDict()

    for rxn_name, reaction in model.reactions.items():
        stoich = [str(int(reaction[species])) for species in model.species_names.values()]
        reaction_set[rxn_name] = f"{{{','.join(stoich)}}}"
    for rxn_name, reaction_reactants in model.reaction_reactants.items():
        reactants_count = [str(int(reaction_reactants[species])) for species in model.species_names.values()]
        reactants_set[rxn_name] = f"{{{','.join(reactants_count)}}}"
    for rxn_name, reaction_products in model.reaction_products.items():
        products_count = [str(int(reaction_products[species])) for species in model.species_names.values()]
        products_set[rxn_name] = f"{{{','.join(products_count)}}}"

    reaction_names = " ".join([f"REACTION_NAME({rxn})" for rxn in reaction_set.keys()])
    reaction_set = f"{{{','.join(reaction_set.values())}}}"
    reactants_set = f"{{{','.join(reactants_set.values())}}}"
    products_set = f"{{{','.join(products_set.values())}}}"

    return {
        "GPY_NUM_REACTIONS": num_reactions,
        "GPY_REACTION_NAMES": reaction_names,
        "GPY_REACTIONS": reaction_set,
        "GPY_REACTION_REACTANTS": reactants_set,
        "GPY_REACTION_PRODUCTS": products_set,
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
