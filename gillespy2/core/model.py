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
import platform
from collections import OrderedDict, ChainMap

import gillespy2
from gillespy2.core.assignmentrule import AssignmentRule
from gillespy2.core.events import Event
from gillespy2.core.functiondefinition import FunctionDefinition
from gillespy2.core.parameter import Parameter
from gillespy2.core.raterule import RateRule
from gillespy2.core.reaction import Reaction
from gillespy2.core.species import Species
from gillespy2.core.timespan import TimeSpan
from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.jsonify import Jsonify, TranslationTable
from gillespy2.core.gillespyError import (
    ModelError,
    SimulationError,
    StochMLImportError,
    InvalidStochMLError,
)

try:
    import lxml.etree as eTree

    no_pretty_print = False

except ImportError:
    import xml.etree.ElementTree as eTree
    import xml.dom.minidom
    import re
    no_pretty_print = True


def import_SBML(filename, name=None, gillespy_model=None, report_silently_with_sbml_error=False):
    """
    SBML to GillesPy model converter. NOTE: non-mass-action rates
    in terms of concentrations may not be converted for population
    simulation. Use caution when importing SBML.

    :param filename: Path to the SBML file for conversion.
    :type filename: str

    :param name: Name of the resulting model
    :type name: str

    :param gillespy_model: If desired, the SBML model may be added to an existing GillesPy model
    :type gillespy_model: gillespy.Model

    :param report_silently_with_sbml_error: SBML import will fail silently and
                SBML errors will not output to the user.
    :type report_silently_with_sbml_error: bool
    """

    try:
        from gillespy2.sbml.SBMLimport import convert # pylint: disable=import-outside-toplevel
    except ImportError as err:
        raise ImportError('SBML conversion not imported successfully') from err

    return convert(
        filename, model_name=name, gillespy_model=gillespy_model,
        report_silently_with_sbml_error=report_silently_with_sbml_error
    )


def export_SBML(gillespy_model, filename=None):
    """
    GillesPy model to SBML converter

    :param gillespy_model: GillesPy model to be converted to SBML
    :type gillespy_model: gillespy.Model

    :param filename: Path to the SBML file for conversion
    :type filename: str
    """
    try:
        from gillespy2.sbml.SBMLexport import export # pylint: disable=import-outside-toplevel
    except ImportError as err:
        raise ImportError('SBML export conversion not imported successfully') from err

    return export(gillespy_model, path=filename)


def export_StochSS(gillespy_model, filename=None, return_stochss_model=False):
    """
    GillesPy model to StochSS converter

    :param gillespy_model: GillesPy model to be converted to StochSS
    :type gillespy_model: gillespy.Model

    :param filename: Path to the StochSS file for conversion
    :type filename: str
    """
    try:
        from gillespy2.stochss.StochSSexport import export # pylint: disable=import-outside-toplevel
    except ImportError as err:
        raise ImportError('StochSS export conversion not imported successfully') from err

    return export(gillespy_model, path=filename, return_stochss_model=return_stochss_model)


class Model(SortableObject, Jsonify):
    """
    Representation of a well mixed biochemical model. Contains reactions,
    parameters, species.

    :param name: The name of the model, or an annotation describing it.
    :type name: str

    :param population: The type of model being described. A discrete stochastic model is a
        population model (True), a deterministic model is a concentration model
        (False). Automatic conversion from population to concentration models
        may be used, by setting the volume parameter.
    :type population: bool

    :param volume: The volume of the system matters when converting to from population to
        concentration form. This will also set a parameter "vol" for use in
        custom (i.e. non-mass-action) propensity functions.
    :type volume: float

    :param tspan: The timepoints at which the model should be simulated. If None, a
        default timespan is added. May be set later, see Model.timespan
    :type tspan: numpy ndarray

    :param annotation: Option further description of model
    :type annotation: str
    """

    # reserved names for model species/parameter names, volume, and operators.
    reserved_names = ['vol']
    special_characters = ['[', ']', '+', '-', '*', '/', '.', '^']

    def __init__(self, name="", population=True, volume=1.0, tspan=None, annotation="model"):
        """ Create an empty model. """

        # The name that the model is referenced by (should be a String)
        self.name = name
        self.annotation = annotation

        # Dictionaries with model element objects.
        # Model element names are used as keys.
        self.listOfParameters = OrderedDict()
        self.listOfSpecies = OrderedDict()
        self.listOfReactions = OrderedDict()

        self.listOfAssignmentRules = OrderedDict()
        self.listOfRateRules = OrderedDict()
        self.listOfEvents = OrderedDict()
        self.listOfFunctionDefinitions = OrderedDict()

        # Dictionaries with model element objects.
        # Model element names are used as keys, and values are
        # sanitized versions of the names/formulas.
        # These dictionaries contain sanitized values and are for
        # Internal use only
        self._listOfParameters = OrderedDict()
        self._listOfSpecies = OrderedDict()
        self._listOfReactions = OrderedDict()
        self._listOfAssignmentRules = OrderedDict()
        self._listOfRateRules = OrderedDict()
        self._listOfEvents = OrderedDict()
        self._listOfFunctionDefinitions = OrderedDict()
        # This defines the unit system at work for all numbers in the model
        # It should be a logical error to leave this undefined, subclasses
        # should set it
        if population:
            self.units = "population"
        else:
            self.units = "concentration"
            if volume != 1.0:
                raise ModelError(
                    "Concentration models account for volume implicitly, explicit volume definition is not required. "
                    "Note: concentration models may only be simulated deterministically."
                )

        self.volume = volume

        # Dict that holds flattended parameters and species for
        # evaluation of expressions in the scope of the model.
        self.namespace = OrderedDict([])

        if tspan is None:
            self.tspan = None
        else:
            self.timespan(tspan)

        # Change Jsonify settings to disable private variable
        # JSON hashing and enable automatic translation table gen.
        self._hash_private_vars = False
        self._generate_translation_table = True

    def __str__(self):
        divider = '\n**********\n'

        def decorate(header):
            return '\n' + divider + header + divider

        print_string = self.name
        if len(self.listOfSpecies):
            print_string += decorate('Species')
            for species in sorted(self.listOfSpecies.values()):
                print_string += '\n' + str(species)
        if len(self.listOfParameters):
            print_string += decorate('Parameters')
            for parameter in sorted(self.listOfParameters.values()):
                print_string += '\n' + str(parameter)
        if len(self.listOfReactions):
            print_string += decorate('Reactions')
            for reaction in sorted(self.listOfReactions.values()):
                print_string += '\n' + str(reaction)
        if len(self.listOfEvents):
            print_string += decorate('Events')
            for event in sorted(self.listOfEvents.values()):
                print_string += '\n' + str(event)
        if len(self.listOfAssignmentRules):
            print_string += decorate('Assignment Rules')
            for assign_rule in sorted(self.listOfAssignmentRules.values()):
                print_string += '\n' + str(assign_rule)
        if len(self.listOfRateRules):
            print_string += decorate('Rate Rules')
            for rate_rule in sorted(self.listOfRateRules.values()):
                print_string += '\n' + str(rate_rule)
        if len(self.listOfFunctionDefinitions):
            print_string += decorate('Function Definitions')
            for func_def in sorted(self.listOfFunctionDefinitions.values()):
                print_string += '\n' + str(func_def)
        if self.tspan is not None:
            print_string += decorate('Timespan')
            print_string += str(self.tspan)
        return print_string

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.get_element(key)
        if hasattr(self.__class__, "__missing__"):
            return self.__class__.__missing__(self, key)
        raise KeyError(f"{key} is an invalid key.")

    def _problem_with_name(self, name):
        if name in Model.reserved_names:
            names = Model.reserved_names
            raise ModelError(
                f'Name "{name}" is unavailable. It is reserved for internal GillesPy use. Reserved Names: ({names}).'
            )
        if name in self.listOfSpecies:
            raise ModelError(f'Name "{name}" is unavailable. A species with that name exists.')
        if name in self.listOfParameters:
            raise ModelError(f'Name "{name}" is unavailable. A parameter with that name exists.')
        if name in self.listOfReactions:
            raise ModelError(f'Name "{name}" is unavailable. A reaction with that name exists.')
        if name in self.listOfEvents:
            raise ModelError(f'Name "{name}" is unavailable. An event with that name exists.')
        if name in self.listOfRateRules:
            raise ModelError(f'Name "{name}" is unavailable. A rate rule with that name exists.')
        if name in self.listOfAssignmentRules:
            raise ModelError(f'Name "{name}" is unavailable. An assignment rule with that name exists.')
        if name in self.listOfFunctionDefinitions:
            raise ModelError(f'Name "{name}" is unavailable. A function definition with that name exists.')
        if name.isdigit():
            raise ModelError(f'Name "{name}" is unavailable. Names must not be numeric strings.')
        for special_character in Model.special_characters:
            if special_character in name:
                chars = Model.special_characters
                raise ModelError(
                    f'Name "{name}" is unavailable. Names must not contain special characters: {chars}.'
                )

    def _resolve_event(self, event):
        def validate(event):
            from gillespy2.core.gillespyError import EventError # pylint: disable=import-outside-toplevel
            if event.trigger is None or not hasattr(event.trigger, 'expression'):
                raise EventError('An Event must contain a valid trigger.')
        try:
            validate(event)

            # Confirm that the variable in the event assignments are part of the model.
            for assign in event.assignments:
                name = assign.variable if isinstance(assign.variable, str) else assign.variable.name
                assign.variable = self.get_element(name)
        except ModelError as err:
            raise ModelError(
                f"Could not add/resolve event: {event.name}, Reason given: {err}"
            ) from err

    def _resolve_all_events(self):
        for _, event in self.listOfEvents.items():
            self._resolve_event(event)

    def _resolve_parameter(self, parameter):
        try:
            parameter.validate()

            # Calculate the parameters value.
            self.update_namespace()
            parameter._evaluate(self.namespace)
        except ModelError as err:
            raise ModelError(f"Could not add/resolve parameter: {parameter.name}, Reason given: {err}") from err

    def _resolve_all_parameters(self):
        for _, parameter in self.listOfParameters.items():
            self._resolve_parameter(parameter)

    def _resolve_rule(self, rule):
        def validate(rule):
            from gillespy2.core.gillespyError import RateRuleError, AssignmentRuleError # pylint: disable=import-outside-toplevel
            errors = {"RateRule": RateRuleError, "AssignmentRule": AssignmentRuleError}
            error_class = errors[type(rule).__name__]
            if rule.variable is None:
                raise error_class('A GillesPy2 Rate/Assignment Rule must be associated with a valid variable')
            if rule.formula == '':
                raise error_class('Invalid Rate/Assignment Rule. Expression must be a non-empty string value')
        try:
            validate(rule)

            # Confirm that the variable is part of the model.
            name = rule.variable if isinstance(rule.variable, str) else rule.variable.name
            rule.variable = self.get_element(name)
        except ModelError as err:
            raise ModelError(
                f"Could not add/resolve rate_rule: {rule.name}, Reason given: {err}"
            ) from err

    def _resolve_all_rate_rules(self):
        for _, rate_rule in self.listOfRateRules.items():
            self._resolve_rule(rate_rule)

    def _resolve_all_assignment_rules(self):
        for _, assign_rule in self.listOfAssignmentRules.items():
            self._resolve_rule(assign_rule)

    def _resolve_reaction(self, reaction):
        try:
            reaction.validate()

            # If the rate parameter exists in the reaction, confirm that it is a part of the model
            if reaction.marate is not None:
                name = reaction.marate if isinstance(reaction.marate, str) else reaction.marate.name
                reaction.marate = self.get_parameter(name)

            # Confirm that all species in reactants are part of the model
            for species in list(reaction.reactants.keys()):
                stoichiometry = reaction.reactants[species]
                name = species if isinstance(species, str) else species.name
                stoich_spec = self.get_species(name)
                if stoich_spec not in reaction.reactants:
                    reaction.reactants[stoich_spec] = stoichiometry
                    del reaction.reactants[species]

            # Confirm that all species in products are part of the model
            for species in list(reaction.products.keys()):
                stoichiometry = reaction.products[species]
                name = species if isinstance(species, str) else species.name
                stoich_spec = self.get_species(name)
                if stoich_spec not in reaction.products:
                    reaction.products[stoich_spec] = stoichiometry
                    del reaction.products[species]
        except ModelError as err:
            raise ModelError(f"Could not add/resolve reaction: {reaction.name}, Reason given: {err}") from err

    def _resolve_all_reactions(self):
        for _, reaction in self.listOfReactions.items():
            self._resolve_reaction(reaction)

    def update_namespace(self):
        """ Create a dict with flattened parameter and species objects. """
        self.namespace = OrderedDict([])
        for param in self.listOfParameters:
            self.namespace[param] = self.listOfParameters[param].value

    def add(self, components):
        """
        Adds a component, or list of components to the model. If a list is provided, Species
        and Parameters are added before other components.  Lists may contain any combination
        of accepted types other than lists and do not need to be in any particular order.

        :param components: The component or list of components to be added the the model.
        :type components: Species, Parameters, Reactions, Events, Rate Rules, Assignment Rules, \
                          FunctionDefinitions, TimeSpan, or list

        :returns: The components that were added to the model.
        :rtype: Species, Parameters, Reactions, Events, Rate Rules, Assignment Rules, \
                FunctionDefinitions, TimeSpan, or list

        :raises ModelError: Component is invalid.
        """
        if isinstance(components, list):
            p_types = (Species, Parameter, FunctionDefinition, TimeSpan)
            p_names = (p_type.__name__ for p_type in p_types)

            others = []
            for component in components:
                if isinstance(component, p_types) or type(component).__name__ in p_names:
                    self.add(component)
                else:
                    others.append(component)

            for component in others:
                self.add(component)
        elif isinstance(components, AssignmentRule) or type(components).__name__ == AssignmentRule.__name__:
            self.add_assignment_rule(components)
        elif isinstance(components, Event) or type(components).__name__ == Event.__name__:
            self.add_event(components)
        elif isinstance(components, FunctionDefinition) or type(components).__name__ == FunctionDefinition.__name__:
            self.add_function_definition(components)
        elif isinstance(components, Parameter) or type(components).__name__ == Parameter.__name__:
            self.add_parameter(components)
        elif isinstance(components, RateRule) or type(components).__name__ == RateRule.__name__:
            self.add_rate_rule(components)
        elif isinstance(components, Reaction) or type(components).__name__ == Reaction.__name__:
            self.add_reaction(components)
        elif isinstance(components, Species) or type(components).__name__ == Species.__name__:
            self.add_species(components)
        elif isinstance(components, TimeSpan) or type(components).__name__ == TimeSpan.__name__:
            self.timespan(components)
        else:
            raise ModelError(f"Unsupported component: {type(components)} is not a valid component.")
        return components

    def get_element(self, name):
        """
        Get a model element specified by name.

        :param name: Name of the element to be returned.
        :type name: str

        :returns: The specified gillespy2.Model element.
        :rtype: Species, Parameters, Reactions, Events, RateRules, AssignmentRules, \
                FunctionDefinitions, or TimeSpan

        :raises ModelError: If the element is not part of the model.
        """
        if name in ("tspan", "timespan"):
            return self.tspan
        if name in self.listOfSpecies:
            return self.get_species(name)
        if name in self.listOfParameters:
            return self.get_parameter(name)
        if name in self.listOfReactions:
            return self.get_reaction(name)
        if name in self.listOfRateRules:
            return self.get_rate_rule(name)
        if name in self.listOfAssignmentRules:
            return self.get_assignment_rule(name)
        if name in self.listOfEvents:
            return self.get_event(name)
        if name in self.listOfFunctionDefinitions:
            return self.get_function_definition(name)
        raise ModelError(f"{self.name} does not contain an element named {name}.")

    def add_species(self, species):
        """
        Adds a species, or list of species to the model.

        :param species: The species or list of species to be added to the model object
        :type species: gillespy2.Species | list of gillespy2.Species

        :returns: The species or list of species that were added to the model.
        :rtype: gillespy2.Species | list of gillespy2.Species

        :raises ModelError: If an invalid species is provided or if Species.validate fails.
        """
        if isinstance(species, list):
            for spec in sorted(species):
                self.add_species(spec)
        elif isinstance(species, Species) or type(species).__name__ == "Species":
            try:
                species.validate()
                self._problem_with_name(species.name)
                self.listOfSpecies[species.name] = species
                self._listOfSpecies[species.name] = f'S{len(self._listOfSpecies)}'
            except ModelError as err:
                errmsg = f"Could not add species: {species.name}, Reason given: {err}"
                raise ModelError(errmsg) from err
        else:
            errmsg = f"species must be of type Species or list of Species not {type(species)}."
            raise ModelError(errmsg)
        return species

    def delete_species(self, name):
        """
        Removes a species object by name.

        :param name: Name of the species object to be removed.
        :type name: str

        :raises ModelError: If the species is not part of the model.
        """
        try:
            self.listOfSpecies.pop(name)
            if name in self._listOfSpecies:
                self._listOfSpecies.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain a species named {name}."
            ) from err

    def delete_all_species(self):
        """
        Removes all species from the model object.
        """
        self.listOfSpecies.clear()
        self._listOfSpecies.clear()

    def get_species(self, name):
        """
        Returns a species object by name.

        :param name: Name of the species object to be returned.
        :type name: str

        :returns: The specified species object.
        :rtype: gillespy2.Species

        :raises ModelError: If the species is not part of the model.
        """
        if name not in self.listOfSpecies:
            raise ModelError(f"{self.name} does not contain a species named {name}.")
        return self.listOfSpecies[name]

    def get_all_species(self):
        """
        Get all of the species in the model object.

        :returns: A dict of all species in the model, in the form: {name : species object}.
        :rtype: OrderedDict
        """
        return self.listOfSpecies

    def sanitized_species_names(self):
        """
        Generate a dictionary mapping user chosen species names to simplified formats which will be used
        later on by GillesPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user species names to their internal GillesPy notation.
        """
        species_name_mapping = OrderedDict([])
        for i, name in enumerate(self.listOfSpecies.keys()):
            species_name_mapping[name] = f'S[{i}]'
        return species_name_mapping

    def add_parameter(self, parameters):
        """
        Adds a parameter, or list of parameters to the model.

        :param parameters:  The parameter or list of parameters to be added to the model object.
        :type parameters: gillespy2.Parameter | list of gillespy2.Parameter

        :returns: A parameter or list of Parameters that were added to the model.
        :rtype: gillespy2.Parameter | list of gillespy2.Parameter

        :raises ModelError: If an invalid parameter is provided or if Parameter.validate fails.
        """
        if isinstance(parameters, list):
            for param in sorted(parameters):
                self.add_parameter(param)
        elif isinstance(parameters, Parameter) or type(parameters).__name__ == 'Parameter':
            self._problem_with_name(parameters.name)
            self._resolve_parameter(parameters)
            self.listOfParameters[parameters.name] = parameters
            self._listOfParameters[parameters.name] = f'P{len(self._listOfParameters)}'
        else:
            errmsg = f"parameters must be of type Parameter or list of Parameter not {type(parameters)}."
            raise ModelError(errmsg)
        return parameters

    def delete_parameter(self, name):
        """
        Removes a parameter object by name.

        :param name: Name of the parameter object to be removed.
        :type name: str

        :raises ModelError: If the parameter is not part of the model.
        """
        try:
            self.listOfParameters.pop(name)
            if name in self._listOfParameters:
                self._listOfParameters.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain a parameter named {name}."
            ) from err

    def delete_all_parameters(self):
        """
        Removes all parameters from model object.
        """
        self.listOfParameters.clear()
        self._listOfParameters.clear()

    def get_parameter(self, name):
        """
        Returns a parameter object by name.

        :param name: Name of the parameter object to be returned
        :type name: str

        :returns: The specified parameter object.
        :rtype: gillespy2.Parameter

        :raises ModelError: If the parameter is not part of the model.
        """
        if name not in self.listOfParameters:
            raise ModelError(f"{self.name} does not contain a parameter named {name}.")
        return self.listOfParameters[name]

    def get_all_parameters(self):
        """
        Get all of the parameters in the model object.

        :returns: A dict of all parameters in the model, in the form: {name : parameter object}
        :rtype: OrderedDict
        """
        return self.listOfParameters

    def sanitized_parameter_names(self):
        """
        Generate a dictionary mapping user chosen parameter names to simplified formats which will be used
        later on by GillesPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user parameter names to their internal GillesPy notation.
        """
        parameter_name_mapping = OrderedDict()
        parameter_name_mapping['vol'] = 'V'
        for i, name in enumerate(self.listOfParameters.keys()):
            if name not in parameter_name_mapping:
                parameter_name_mapping[name] = f'P{i}'
        return parameter_name_mapping

    def add_reaction(self, reactions):
        """
        Adds a reaction, or list of reactions to the model.

        :param reactions: The reaction or list of reactions to be added to the model object
        :type reactions: gillespy2.Reaction | list of gillespy2.Reaction

        :returns: The reaction or list of reactions that were added to the model.
        :rtype: gillespy2.Reaction | list of gillespy2.Reaction

        :raises ModelError: If an invalid reaction is provided or if Reaction.validate fails.
        """
        if isinstance(reactions, list):
            for reaction in sorted(reactions):
                self.add_reaction(reaction)
        elif isinstance(reactions, Reaction) or type(reactions).__name__ == "Reaction":
            self._problem_with_name(reactions.name)
            self._resolve_reaction(reactions)
            self.listOfReactions[reactions.name] = reactions
            # Build Sanitized reaction as well
            sanitized_reaction = reactions._create_sanitized_reaction(
                len(self.listOfReactions), self._listOfSpecies, self._listOfParameters
            )
            self._listOfReactions[reactions.name] = sanitized_reaction
        else:
            errmsg = f"reactions must be of type Reaction or list of Reaction not {type(reactions)}."
            raise ModelError(errmsg)
        return reactions

    def delete_reaction(self, name):
        """
        Removes a reaction object by name.

        :param name: Name of the reaction object to be removed.
        :type name: str

        :raises ModelError: If the reaction is not part of the model.
        """
        try:
            self.listOfReactions.pop(name)
            if name in self._listOfReactions:
                self._listOfReactions.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain a reaction named {name}."
            ) from err

    def delete_all_reactions(self):
        """
        Removes all reactions from the model object.
        """
        self.listOfReactions.clear()
        self._listOfReactions.clear()

    def get_reaction(self, name):
        """
        Returns a reaction object by name.

        :param name: Name of the reaction object to be returned
        :type name: str

        :returns: The specified reaction object.
        :rtype: gillespy2.Reaction

        :raises ModelError: If the reaction is not part of the model.
        """
        if name not in self.listOfReactions:
            raise ModelError(f"{self.name} does not contain a reaction named {name}.")
        return self.listOfReactions[name]

    def get_all_reactions(self):
        """
        Get all of the reactions in the model object.

        :returns: A dict of all reactions in the model, in the form: {name : reaction object}.
        :rtype: OrderedDict
        """
        return self.listOfReactions

    def add_rate_rule(self, rate_rule):
        """
        Adds a rate rule, or list of rate rules to the model.

        :param rate_rule: The rate rule or list of rate rules to be added to the model object.
        :type rate_rule: gillespy2.RateRule | list of gillespy2.RateRules

        :returns: The rate rule or list of rate rules that were added to the model.
        :rtype: gillespy2.RateRule | list of gillespy2.RateRule

        :raises ModelError: If an invalid rate rule is provided or if rate rule validation fails.
        """
        if isinstance(rate_rule, list):
            for r_rule in sorted(rate_rule):
                self.add_rate_rule(r_rule)
        elif isinstance(rate_rule, RateRule) or type(rate_rule).__name__ == "RateRule":
            self._problem_with_name(rate_rule.name)
            ar_vars = [a_rule.variable for a_rule in self.listOfAssignmentRules.values()]
            rr_vars = [r_rule.variable for r_rule in self.listOfRateRules.values()]
            if rate_rule.variable in ar_vars:
                raise ModelError(
                    f"Duplicate variable in rate_rules AND assignment_rules: {rate_rule.variable}."
                )
            if rate_rule.variable in rr_vars:
                raise ModelError(f"Duplicate variable in rate_rules: {rate_rule.variable}.")
            self._resolve_rule(rate_rule)
            if rate_rule.variable.name in self.listOfSpecies:
                # check if the rate_rule's target's mode is continious
                if rate_rule.variable.mode == 'discrete':
                    raise ModelError("RateRules can not target discrete species")
                if rate_rule.variable.mode is None or rate_rule.variable.mode == 'dynamic':
                    rate_rule.variable.mode = 'continuous'
                    from gillespy2.core import log # pylint:disable=import-outside-toplevel
                    errmsg = f"Changing {rate_rule.variable.name}.mode='continuous' as it is the target of RateRule"
                    errmsg += f" {rate_rule.name}"
                    log.warning(errmsg)

            self.listOfRateRules[rate_rule.name] = rate_rule
            # Build the sanitized rate rule
            sanitized_rate_rule = RateRule(name=f'RR{len(self._listOfRateRules)}')
            sanitized_rate_rule.formula = rate_rule.sanitized_formula(
                self._listOfSpecies, self._listOfParameters
            )
            self._listOfRateRules[rate_rule.name] = sanitized_rate_rule
        else:
            errmsg = f"rate_rule must be of type RateRule or list of RateRules not {type(rate_rule)}."
            raise ModelError(errmsg)
        return rate_rule

    def delete_rate_rule(self, name):
        """
        Removes rate rule object by name.

        :param name: Name of the rate rule to be removed.
        :type name: str

        :raises ModelError: If the rate rule is not part of the model.
        """
        try:
            self.listOfRateRules.pop(name)
            if name in self._listOfRateRules:
                self._listOfRateRules.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain a rate rule named {name}."
            ) from err

    def delete_all_rate_rules(self):
        """
        Removes all rate rules from the model object.
        """
        self.listOfRateRules.clear()
        self._listOfRateRules.clear()

    def get_rate_rule(self, name):
        """
        Returns a rate rule object by name.

        :param name: Name of the rate rule object to be returned.
        :type name: str

        :returns: The specified rate rule object.
        :rtype: gillespy2.RateRule

        :raises ModelError: If the rate rule is not part of the model.
        """
        if name not in self.listOfRateRules:
            raise ModelError(f"{self.name} does not contain a rate rule named {name}.")
        return self.listOfRateRules[name]

    def get_all_rate_rules(self):
        """
        Get all of the rate rules in the model object.

        :returns: A dict of all rate rules in the model, in the form: {name : rate rule object}.
        :rtype: OrderedDict
        """
        return self.listOfRateRules

    def add_assignment_rule(self, assignment_rule):
        """
        Add an assignment rule, or list of assignment rules to the model.

        :param assignment_rules: The assignment rule or list of assignment rules to be added to the model object.
        :type assignment_rules: gillespy2.AssignmentRule or list of gillespy2.AssignmentRules

        :returns: The assignment rule or list of assignment rules that were added to the model.
        :rtype: gillespy2.AssignmentRule | list of gillespy2.AssignmentRule

        :raises ModelError: If an invalid assignment rule is provided or if assignment rule validation fails.
        """
        if isinstance(assignment_rule, list):
            for a_rule in assignment_rule:
                self.add_assignment_rule(a_rule)
        elif isinstance(assignment_rule, AssignmentRule) or type(assignment_rule).__name__ == "AssignmentRule":
            self._problem_with_name(assignment_rule.name)
            ar_vars = [a_rule.variable for a_rule in self.listOfAssignmentRules.values()]
            rr_vars = [r_rule.variable for r_rule in self.listOfRateRules.values()]
            if assignment_rule.variable in rr_vars:
                raise ModelError(
                    f"Duplicate variable in rate_rules AND assignment_rules: {assignment_rule.variable}."
                )
            if assignment_rule.variable in ar_vars:
                raise ModelError(f"Duplicate variable in assignments_rules: {assignment_rule.variable}.")
            self._resolve_rule(assignment_rule)
            self.listOfAssignmentRules[assignment_rule.name] = assignment_rule
            # Build the sanitized assignment rule
            sanitized_assignment_rule = AssignmentRule(name=f'AR{len(self._listOfAssignmentRules)}')
            sanitized_assignment_rule.formula = assignment_rule.sanitized_formula(
                self._listOfSpecies, self._listOfParameters
            )
            self._listOfAssignmentRules[assignment_rule.name] = sanitized_assignment_rule
        else:
            errmsg = "assignment_rule must be of type AssignmentRule or "
            errmsg += f"list of AssignmentRules not {type(assignment_rule)}."
            raise ModelError(errmsg)
        return assignment_rule

    def delete_assignment_rule(self, name):
        """
        Removes an assignment rule object by model.

        :param name: Name of the assignment rule object to be removed.
        :type name: str

        :raises ModelError: If the assignment rule is not part of the model.
        """
        try:
            self.listOfAssignmentRules.pop(name)
            if name in self._listOfAssignmentRules:
                self._listOfAssignmentRules.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain an assignment rule named {name}."
            ) from err

    def delete_all_assignment_rules(self):
        """
        Removes all assignment rules from the model object.
        """
        self.listOfAssignmentRules.clear()
        self._listOfAssignmentRules.clear()

    def get_assignment_rule(self, name):
        """
        Returns an assignment rule object by name.

        :param name: Name of the assignment rule object to be returned.
        :type name: str

        :returns: The specified assignment rule object.
        :rtype: gillespy2.AssignmentRule

        :raises ModelError: If the assignment rule is not part of the model.
        """
        if name not in self.listOfAssignmentRules:
            raise ModelError(f"{self.name} does not contain an assignment rule named {name}.")
        return self.listOfAssignmentRules[name]

    def get_all_assignment_rules(self):
        """
        Get all of the assignment rules in the model object.

        :returns: A dict of all assignemt rules in the model, in the form: {name: reaction object}.
        """
        return self.listOfAssignmentRules

    def add_event(self, event):
        """
        Adds an event, or list of events to the model.

        :param event: The event or list of event to be added to the model object.
        :type event: gillespy2.Event | list of gillespy2.Events

        :returns: The event or list of events that were added to the model.
        :rtype: gillespy2.Event | list of gillespy2.Event

        :raises ModelError: If an invalid event is provided or if event validation fails.
        """
        if isinstance(event, list):
            for evnt in event:
                self.add_event(evnt)
        elif isinstance(event, Event) or type(event).__name__ == "Event":
            self._problem_with_name(event.name)
            self._resolve_event(event)
            self.listOfEvents[event.name] = event
        else:
            errmsg = f"event must be of type Event or list of Events not {type(event)}"
            raise ModelError(errmsg)
        return event

    def delete_event(self, name):
        """
        Removes an event object by name.

        :param name: Name of the event object to be removed.
        :type name: str
        """
        try:
            self.listOfEvents.pop(name)
            if name in self._listOfEvents:
                self._listOfEvents.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain an event named {name}."
            ) from err

    def delete_all_events(self):
        """
        Removes all events from the model object.
        """
        self.listOfEvents.clear()
        self._listOfEvents.clear()

    def get_event(self, name):
        """
        Returns an event object by name.

        :param name: Name of the event object to be returned.
        :type name: str

        :returns: The specified event object.
        :rtype: gillespy2.Event

        :raises ModelError: If the event is not part of the model.
        """
        if name not in self.listOfEvents:
            raise ModelError(f"{self.name} does not contain an event named {name}.")
        return self.listOfEvents[name]

    def get_all_events(self):
        """
        Get all of the events in the model object.

        :returns: A dict of all evetns in the model, in the form: {name : event object}
        """
        return self.listOfEvents

    def add_function_definition(self, function_definition):
        """
        Add function definition, or list of function definitions to the model

        :param function_definition: The function definition, or list of function definitions \
                to be added to the model object.
        :type function_definition: gillespy2.FunctionDefinition | list of gillespy2.FunctionDefinitions.

        :returns: The function defintion or list of function definitions that were added to the model.
        :rtype: gillespy2.FunctionDefinitions | list of gillespy2.FunctionDefinitions

        :raises ModelError: If an invalid function definition is provided.
        """
        if isinstance(function_definition, list):
            for func_def in function_definition:
                self.add_function_definition(func_def)
        elif isinstance(function_definition, FunctionDefinition) or \
                        type(function_definition).__name__ == "FunctionDefinition":
            self._problem_with_name(function_definition.name)
            self.listOfFunctionDefinitions[function_definition.name] = function_definition
        else:
            errmsg = "function_definition must be of type FunctionDefinition or "
            errmsg += f"list of FunctionDefinitions not {type(function_definition)}."
            raise ModelError(errmsg)

    def delete_function_definition(self, name):
        """
        Removes a function definition object by name.

        :param name: Name of the function definition object to be removed.
        :type name: str
        """
        try:
            self.listOfFunctionDefinitions.pop(name)
            if name in self._listOfFunctionDefinitions:
                self._listOfFunctionDefinitions.pop(name)
        except KeyError as err:
            raise ModelError(
                f"{self.name} does not contain a function definition named {name}."
            ) from err

    def delete_all_function_definitions(self):
        """
        Removes all function definitions from the model object.
        """
        self.listOfFunctionDefinitions.clear()
        self._listOfFunctionDefinitions.clear()

    def get_function_definition(self, name):
        """
        Returns a function definition object by name.

        :param name: Name of the function definition object to be returned.
        :type name: str

        :returns: The specified function definition object.
        :rtype: gillespy2.FunctionDefinition
        """
        if name not in self.listOfFunctionDefinitions:
            raise ModelError(f"{self.name} does not contain a function definition named {name}.")
        return self.listOfFunctionDefinitions[name]

    def get_all_function_definitions(self):
        """
        Get all of the function definitions in the model object.

        :returns: A dict of all function definitions in the model, in the form {name : function definition object}.
        :rtype: OrderedDict
        """
        return self.listOfFunctionDefinitions

    def timespan(self, time_span):
        """
        Set the time span of simulation. StochKit does not support non-uniform
        timespans.

        :param time_span: Evenly-spaced list of times at which to sample the species populations
            during the simulation. Best to use the form gillespy2.TimeSpan(np.linspace(
            <start time>, <end time>, <number of time-points, inclusive>))
        :type time_span: gillespy2.TimeSpan | iterator
        """
        if isinstance(time_span, TimeSpan) or type(time_span).__name__ == "TimeSpan":
            self.tspan = time_span
        else:
            self.tspan = TimeSpan(time_span)

    def make_translation_table(self):
        species = self.listOfSpecies.values()
        reactions = self.listOfReactions.values()
        parameters = self.listOfParameters.values()
        assignments = self.listOfAssignmentRules.values()
        rates = self.listOfRateRules.values()
        events = self.listOfEvents.values()
        for event in events:
            for assignment in event.__dict__['assignments']:
                del assignment.__dict__['_EventAssignment__name_deprecated']
        functions = self.listOfFunctionDefinitions.values()

        # A translation table is used to anonymize user-defined variable names and formulas into generic counterparts.
        translation_table = dict(ChainMap(

            # Build translation mappings for user-defined variable names.
            dict({ self.name: "Model" }),
            dict(zip((str(x.name) for x in species), (f"S_{x + 100}" for x in range(0, len(species))))),
            dict(zip((str(x.name) for x in reactions), (f"R_{x + 100}" for x in range(0, len(reactions))))),
            dict(zip((str(x.name) for x in parameters), (f"P_{x + 100}" for x in range(0, len(parameters))))),
            dict(zip((str(x.name) for x in assignments), (f"AR_{x + 100}" for x in range(0, len(assignments))))),
            dict(zip((str(x.name) for x in rates), (f"RR_{x + 100}" for x in range(0, len(rates))))),
            dict(zip((str(x.name) for x in events), (f"E_{x + 100}" for x in range(0, len(events))))),
            dict(zip((str(x.name) for x in functions), (f"F_{x + 100}" for x in range(0, len(functions))))),
        ))

        return TranslationTable(to_anon=translation_table)

    def serialize(self):
        """ Serializes the Model object to valid StochML. """
        doc = StochMLDocument().from_model(self)
        return doc.to_string()

    def set_units(self, units):
        """
        Sets the units of the model to either "population" or "concentration"

        :param units: Either "population" or "concentration"
        :type units: str
        """
        if units.lower() == 'concentration' or units.lower() == 'population':
            self.units = units.lower()
        else:
            raise ModelError("units must be either concentration or population (case insensitive)")

    def get_best_solver(self):
        """
        Finds best solver for the users simulation. Currently, AssignmentRules, RateRules, FunctionDefinitions,
        Events, and Species with a dynamic, or continuous population must use the TauHybridSolver.

        :param precompile: If True, and the model contains no AssignmentRules, RateRules, FunctionDefinitions, Events,
            or Species with a dynamic or continuous population, it will choose SSACSolver
        :type precompile: bool

        :returns: gillespy2.gillespySolver
        """
        hybrid_check = False
        if len(self.get_all_rate_rules())  or len(self.get_all_events()):
            hybrid_check = True

        if len(self.get_all_species()) and not hybrid_check:
            for i in self.get_all_species():
                tempMode = self.get_species(i).mode
                if tempMode in ('dynamic', 'continuous'):
                    hybrid_check = True
                    break

        chybrid_check = hybrid_check
        if len(self.get_all_assignment_rules()) or len(self.get_all_function_definitions()):
            hybrid_check = True
            chybrid_check = False

        from gillespy2.solvers.cpp.build.build_engine import BuildEngine # pylint: disable=import-outside-toplevel
        missing_deps = BuildEngine.get_missing_dependencies()
        windows_space = platform.system() == "Windows" and " " in gillespy2.__file__
        can_use_cpp = len(missing_deps) == 0 and not windows_space

        if hybrid_check:
            if can_use_cpp and chybrid_check:
                from gillespy2 import TauHybridCSolver # pylint: disable=import-outside-toplevel
                return TauHybridCSolver

            from gillespy2 import TauHybridSolver # pylint: disable=import-outside-toplevel
            return TauHybridSolver

        if can_use_cpp:
            from gillespy2 import SSACSolver # pylint: disable=import-outside-toplevel
            return SSACSolver

        from gillespy2 import NumPySSASolver # pylint: disable=import-outside-toplevel
        return NumPySSASolver

    def get_best_solver_algo(self, algorithm):
        """
        If user has specified a particular algorithm, we return either the Python or C++ version of that algorithm
        """
        from gillespy2.solvers.cpp.build.build_engine import BuildEngine # pylint: disable=import-outside-toplevel
        missing_deps = BuildEngine.get_missing_dependencies()
        windows_space = platform.system() == "Windows" and " " in gillespy2.__file__
        can_use_cpp = len(missing_deps) == 0 and not windows_space
        chybrid_check = not (len(self.get_all_assignment_rules()) or len(self.get_all_function_definitions()))

        if algorithm == 'Tau-Leaping':
            if can_use_cpp:
                from gillespy2 import TauLeapingCSolver # pylint: disable=import-outside-toplevel
                return TauLeapingCSolver

            from gillespy2 import TauLeapingSolver # pylint: disable=import-outside-toplevel
            return TauLeapingSolver

        if algorithm == 'SSA':
            if can_use_cpp:
                from gillespy2 import SSACSolver # pylint: disable=import-outside-toplevel
                return SSACSolver

            from gillespy2 import NumPySSASolver # pylint: disable=import-outside-toplevel
            return NumPySSASolver

        if algorithm == 'ODE':
            if can_use_cpp:
                from gillespy2 import ODECSolver # pylint: disable=import-outside-toplevel
                return ODECSolver

            from gillespy2 import ODESolver # pylint: disable=import-outside-toplevel
            return ODESolver

        if algorithm == 'Tau-Hybrid':
            if can_use_cpp and chybrid_check:
                from gillespy2 import TauHybridCSolver # pylint: disable=import-outside-toplevel
                return TauHybridCSolver

            from gillespy2 import TauHybridSolver # pylint: disable=import-outside-toplevel
            return TauHybridSolver

        if algorithm == 'CLE':
            from gillespy2 import CLESolver # pylint: disable=import-outside-toplevel
            return CLESolver

        raise ModelError("Invalid value for the argument 'algorithm' entered. "
                         "Please enter 'SSA', 'ODE', 'CLE', 'Tau-Leaping', or 'Tau-Hybrid'.")

    def get_model_features(self) -> "Set[Type]":
        """
        Determine what solver-specific model features are present on the model.
        Used to validate that the model is compatible with the given solver.

        :returns: Set containing the classes of every solver-specific feature present on the model.
        """
        features = set()
        if len(self.listOfEvents):
            features.add(gillespy2.Event)
        if len(self.listOfRateRules):
            features.add(gillespy2.RateRule)
        if len(self.listOfAssignmentRules):
            features.add(gillespy2.AssignmentRule)
        if len(self.listOfFunctionDefinitions):
            features.add(gillespy2.FunctionDefinition)
        return features

    def compile_prep(self):
        """
        Prepare the model for export or simulation.
        """
        for _, species in self.listOfSpecies.items():
            species.validate()
        self._resolve_all_parameters()
        self._resolve_all_reactions()
        self._resolve_all_rate_rules()
        self._resolve_all_assignment_rules()
        self._resolve_all_events()

        if self.tspan is not None:
            if not isinstance(self.tspan, TimeSpan) or type(self.tspan).__name__ != "TimeSpan":
                tspan = TimeSpan(self.tspan)
                self.timespan(tspan)
            else:
                self.tspan.validate()

    def run(self, solver=None, timeout=0, t=None, increment=None, algorithm=None, **solver_args):
        """
        Function calling simulation of the model. There are a number of
        parameters to be set here.

        :param solver: The solver by which to simulate the model. This solver object may
            be initialized separately to specify an algorithm. Optional, defaults to ssa solver.
        :type solver: gillespy.GillesPySolver

        :param timeout: Allows a time_out value in seconds to be sent to a signal handler,
            restricting simulation run-time
        :type timeout: int

        :param t: End time of simulation
        :type t: int

        :param solver_args: Solver-specific arguments to be passed to solver.run()

        :param algorithm: Specify algorithm ('ODE', 'Tau-Leaping', or 'SSA') for GillesPy2 to automatically
            pick best solver using that algorithm.
        :type algorithm: str

        :returns:  Returns a Results object that inherits UserList and contains one or more Trajectory objects that
            inherit UserDict. Results object supports graphing and csv export.

        To pause a simulation and retrieve data before the simulation, keyboard interrupt the simulation by pressing
        control+c or pressing stop on a jupyter notebook. To resume a simulation, pass your previously ran results
        into the run method, and set t = to the time you wish the resuming simulation to end (run(resume=results, t=x)).

        **Pause/Resume is only supported for SINGLE TRAJECTORY simulations.
            T MUST BE SET OR UNEXPECTED BEHAVIOR MAY OCCUR.**
        """

        if solver is None:
            if algorithm is not None:
                solver = self.get_best_solver_algo(algorithm)
            else:
                solver = self.get_best_solver()

        if not hasattr(solver, "is_instantiated"):
            try:
                sol_kwargs = {'model': self}
                if "CSolver" in solver.name and \
                    ("resume" in solver_args or "variables" in solver_args or "live_output" in solver_args):
                    sol_kwargs['variable'] = True
                solver = solver(**sol_kwargs)
            except Exception as err:
                raise SimulationError(f"{solver} is not a valid solver.  Reason Given: {err}.") from err

        try:
            return solver.run(t=t, increment=increment, timeout=timeout, **solver_args)
        except Exception as err:
            raise SimulationError(
                f"argument 'solver={solver}' to run() failed.  Reason Given: {err}"
            ) from err

    def problem_with_name(self, *args):
        """
        (deprecated)
        """
        from gillespy2.core import log # pylint: disable=import-outside-toplevel
        log.warning(
            """
            Model.problem_with_name has been deprecated.  Future releases of GillesPy2 may
            not support this feature.  This is an internal function and should not be used.
            """
        )
        self._problem_with_name(*args)

    def set_parameter(self, p_name, expression):
        """
        Set the value of an existing parameter "pname" to "expression" (deprecated).

        :param p_name: Name of the parameter whose value will be set.
        :type p_name: str

        :param expression: String that may be executed in C, describing the value of the
            parameter. May reference other parameters by name. (e.g. "k1*4")
        :type expression: str
        """
        from gillespy2.core import log # pylint: disable=import-outside-toplevel
        log.warning(
            """
            Model.set_parameter has been deprecated.  Future releases of GillesPy2 may
            not support this feature.  Parameter.expression should only be set in the constructor.
            """
        )

        parameter = self.listOfParameters[p_name]
        parameter.expression = expression
        self._resolve_parameter(parameter)

    def validate_reactants_and_products(self, reactions):
        """
        Internal function (deprecated):
        Ensure that the rate and all reactants and products are present in the model
        for the given reaction.  This methods must be called before exporting the model.

        :param reaction: The target reaction to resolve.
        :type reaction: gillespy2.Reaction

        :raises ModelError: If the reaction can't be resolved.
        """
        from gillespy2.core import log # pylint: disable=import-outside-toplevel
        log.warning(
            """
            Model.validate_reactants_and_products has been deprecated. Future releases of
            GillesPy2 may not support this feature.  Use Model._resolve_reaction instead.
            """
        )

        self._resolve_reaction(reactions)

class StochMLDocument():
    """ Serializiation and deserialization of a Model to/from
        the native StochKit2 XML format. """

    def __init__(self):
        # The root element
        self.document = eTree.Element("Model")
        self.annotation = None

    @classmethod
    def from_model(cls, model):
        """
        Creates an StochKit XML document from an exisiting Model object.
        This method assumes that all the parameters in the model are already
        resolved to scalar floats (see Model.resolveParamters).

        Note, this method is intended to be used internally by the models
        'serialization' function, which performs additional operations and
        tests on the model prior to writing out the XML file.

        You should NOT do:

        .. code-block:: python

            document = StochMLDocument.fromModel(model)
            print document.toString()

        You SHOULD do:

        .. code-block:: python

            print model.serialize()

        """

        # Description
        sml_model = cls()

        descript = eTree.Element('Description')

        # Prepare model for export
        model.compile_prep()

        if model.units.lower() == "concentration":
            descript.set('units', model.units.lower())

        descript.text = model.annotation
        sml_model.document.append(descript)

        # Number of Reactions
        num_reacs = eTree.Element('NumberOfReactions')
        num_reacs.text = str(len(model.listOfReactions))
        sml_model.document.append(num_reacs)

        # Number of Species
        num_specs = eTree.Element('NumberOfSpecies')
        num_specs.text = str(len(model.listOfSpecies))
        sml_model.document.append(num_specs)

        # Species
        spec = eTree.Element('SpeciesList')
        for sname in model.listOfSpecies:
            spec.append(sml_model.__species_to_element(model.listOfSpecies[sname]))
        sml_model.document.append(spec)

        # Parameters
        params = eTree.Element('ParametersList')
        for pname in model.listOfParameters:
            params.append(sml_model.__parameter_to_element(
                model.listOfParameters[pname]))

        vol = Parameter(name='vol', expression=model.volume)
        vol._evaluate()
        params.append(sml_model.__parameter_to_element(vol))

        sml_model.document.append(params)

        # Reactions
        reacs = eTree.Element('ReactionsList')
        for rname in model.listOfReactions:
            reacs.append(sml_model.__reaction_to_element(model.listOfReactions[rname], model.volume))
        sml_model.document.append(reacs)

        return sml_model

    @classmethod
    def from_file(cls, filepath):
        """ Intializes the document from an exisiting native StochKit XML
        file read from disk. """
        tree = eTree.parse(filepath)
        root = tree.getroot()
        sml_model = cls()
        sml_model.document = root
        return sml_model

    @classmethod
    def from_string(cls, string):
        """ Intializes the document from an exisiting native StochKit XML
        file read from disk. """
        root = eTree.fromString(string)

        sml_model = cls()
        sml_model.document = root
        return sml_model

    def to_model(self, name):
        """ Instantiates a Model object from a StochMLDocument. """

        # Empty model
        model = Model(name=name)
        root = self.document

        # Try to set name from document
        if model.name == "":
            name = root.find('Name')
            if name.text is None:
                raise NameError("The Name cannot be none")
            model.name = name.text

        # Set annotiation
        ann = root.find('Description')
        if ann is not None:
            units = ann.get('units')

            if units:
                units = units.strip().lower()

            if units == "concentration":
                model.units = "concentration"
            elif units == "population":
                model.units = "population"
            else:  # Default
                model.units = "population"

            if ann.text is None:
                model.annotation = ""
            else:
                model.annotation = ann.text

        # Set units
        units = root.find('Units')
        if units is not None:
            if units.text.strip().lower() == "concentration":
                model.units = "concentration"
            elif units.text.strip().lower() == "population":
                model.units = "population"
            else:  # Default
                model.units = "population"

        # Create parameters
        for sml_param in root.iter('Parameter'):
            name = sml_param.find('Id').text
            expr = sml_param.find('Expression').text
            if name.lower() == 'vol' or name.lower() == 'volume':
                model.volume = float(expr)
            else:
                param = Parameter(name, expression=expr)
                # Try to evaluate the expression in the empty namespace
                # (if the expr is a scalar value)
                param._evaluate()
                model.add_parameter(param)

        # Create species
        for sml_spec in root.iter('Species'):
            name = sml_spec.find('Id').text
            val = sml_spec.find('InitialPopulation').text
            if '.' in val:
                val = float(val)
            else:
                val = int(val)
            spec = Species(name, initial_value=val)
            model.add_species([spec])

        # The namespace_propensity for evaluating the propensity function
        # for reactions must contain all the species and parameters.
        namespace_propensity = OrderedDict()
        all_species = model.get_all_species()
        all_parameters = model.get_all_parameters()

        for param in all_species:
            namespace_propensity[param] = all_species[param].initial_value

        for param in all_parameters:
            namespace_propensity[param] = all_parameters[param].value

        # Create reactions
        for reac in root.iter('Reaction'):
            try:
                name = reac.find('Id').text
            except Exception as err:
                raise InvalidStochMLError("Reaction has no name.") from err

            # Type may be 'mass-action','customized'
            try:
                r_type = reac.find('Type').text
            except Exception as err:
                raise InvalidStochMLError("No reaction type specified.") from err

            g_reactants = {}
            reactants = reac.find('Reactants')
            try:
                for stoich_spec in reactants.iter('SpeciesReference'):
                    specname = stoich_spec.get('id')
                    # The stochiometry should be an integer value, but some
                    # exising StoxhKit models have them as floats. This is
                    # why we need the slightly odd conversion below.
                    stoch = int(float(stoich_spec.get('stoichiometry')))
                    # Select a reference to species with name specname
                    sref = model.listOfSpecies[specname]
                    try:
                        # The sref list should only contain one element if
                        # the XML file is valid.
                        g_reactants[sref] = stoch
                    except Exception as err:
                        raise StochMLImportError(f"Reason given: {err}") from err
            except Exception:
                # Yes, this is correct. 'reactants' can be None
                pass

            g_products = {}
            products = reac.find('Products')
            try:
                for stoich_spec in products.iter('SpeciesReference'):
                    specname = stoich_spec.get('id')
                    stoch = int(float(stoich_spec.get('stoichiometry')))
                    sref = model.listOfSpecies[specname]
                    try:
                        # The sref list should only contain one element if
                        # the XML file is valid.
                        g_products[sref] = stoch
                    except Exception as err:
                        raise StochMLImportError(f"Reason given: {err}") from err
            except Exception:
                # Yes, this is correct. 'products' can be None
                pass

            kwargs = {}
            if r_type == 'mass-action':
                # If it is mass-action, a parameter reference is needed.
                # This has to be a reference to a species instance. We
                # explicitly disallow a scalar value to be passed as the
                # parameter.
                try:
                    ratename = reac.find('Rate').text
                    try:
                        kwargs['rate'] = model.listOfParameters[ratename]
                    except KeyError:
                        # No paramter name is given. This is a valid use case
                        # in StochKit. We generate a name for the paramter,
                        # and create a new parameter instance. The parameter's
                        # value should now be found in 'ratename'.
                        generated_rate_name = "Reaction_" + name + \
                                              "_rate_constant"
                        param = Parameter(name=generated_rate_name,
                                      expression=ratename)
                        # Try to evaluate the parameter to set its value
                        model.add_parameter(param)
                        kwargs['rate'] = model.listOfParameters[generated_rate_name]
                except Exception as err:
                    raise StochMLImportError(f"Reason given: {err}") from err
            elif r_type == 'customized':
                try:
                    propfunc = reac.find('PropensityFunction').text
                except Exception as err:
                    raise InvalidStochMLError(
                        f"Found a customized propensity function, but no expression was given. {err}"
                    ) from err
                kwargs['propensity_function'] = propfunc
            else:
                raise InvalidStochMLError(
                    "Unsupported or no reaction type given for reaction" + name)

            reaction = Reaction(
                name=name, reactants=g_reactants, products=g_products, **kwargs
            )

            model.add_reaction(reaction)

        return model

    def to_string(self):
        """ Returns  the document as a string. """
        try:
            doc = eTree.tostring(self.document, pretty_print=True)
            return doc.decode("utf-8")
        except Exception:
            # Hack to print pretty xml without pretty-print
            # (requires the lxml module).
            doc = eTree.tostring(self.document)
            xmldoc = xml.dom.minidom.parseString(doc)
            ugly_xml = xmldoc.toprettyxml(indent='  ')
            text_re = re.compile(">\n\s+([^<>\s].*?)\n\s+</", re.DOTALL)
            pretty_xml = text_re.sub(">\g<1></", ugly_xml)
            return pretty_xml

    def __species_to_element(self, species):
        sml_spec = eTree.Element('Species')
        id_element = eTree.Element('Id')
        id_element.text = species.name
        sml_spec.append(id_element)

        if hasattr(species, 'description'):
            description_element = eTree.Element('Description')
            description_element.text = species.description
            sml_spec.append(description_element)

        initial_population_element = eTree.Element('InitialPopulation')
        initial_population_element.text = str(species.initial_value)
        sml_spec.append(initial_population_element)

        return sml_spec

    def __parameter_to_element(self, parameter):
        sml_param = eTree.Element('Parameter')
        id_element = eTree.Element('Id')
        id_element.text = parameter.name
        sml_param.append(id_element)
        expression_element = eTree.Element('Expression')
        expression_element.text = str(parameter.value)
        sml_param.append(expression_element)
        return sml_param

    def __reaction_to_element(self, reaction, model_volume):
        sml_reac = eTree.Element('Reaction')

        id_element = eTree.Element('Id')
        id_element.text = reaction.name
        sml_reac.append(id_element)

        description_element = eTree.Element('Description')
        description_element.text = self.annotation
        sml_reac.append(description_element)

        # StochKit2 wants a rate for mass-action propensites
        if reaction.massaction and model_volume == 1.0:
            rate_element = eTree.Element('Rate')
            # A mass-action reactions should only have one parameter
            rate_element.text = reaction.marate.name
            type_element = eTree.Element('Type')
            type_element.text = 'mass-action'
            sml_reac.append(type_element)
            sml_reac.append(rate_element)

        else:
            type_element = eTree.Element('Type')
            type_element.text = 'customized'
            sml_reac.append(type_element)
            function_element = eTree.Element('PropensityFunction')
            function_element.text = reaction.propensity_function
            sml_reac.append(function_element)

        reactants = eTree.Element('Reactants')

        for reactant, stoichiometry in reaction.reactants.items():
            sr_element = eTree.Element('SpeciesReference')
            sr_element.set('id', str(reactant.name))
            sr_element.set('stoichiometry', str(stoichiometry))
            reactants.append(sr_element)

        sml_reac.append(reactants)

        products = eTree.Element('Products')
        for product, stoichiometry in reaction.products.items():
            sr_element = eTree.Element('SpeciesReference')
            sr_element.set('id', str(product.name))
            sr_element.set('stoichiometry', str(stoichiometry))
            products.append(sr_element)
        sml_reac.append(products)

        return sml_reac
