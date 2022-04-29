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
import numpy as np
from collections import OrderedDict

from typing import (
    List,
    Dict,
    Tuple,
    Union,
    Optional,
    Iterable,
    AbstractSet
)

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
from gillespy2.core.results import Trajectory, Results
from gillespy2.core.gillespySolver import GillesPySolver
from gillespy2.core.gillespyError import (
    ParameterError,
    ModelError,
    SimulationError,
    StochMLImportError,
    InvalidStochMLError
)

try:
    import lxml.etree as eTree

    no_pretty_print = False

except:
    import xml.etree.ElementTree as eTree
    import xml.dom.minidom
    import re
    no_pretty_print = True

# A component can be one or more of the following types.
Component = Union[
    Event,
    Species,
    TimeSpan,
    Reaction,
    RateRule,
    Parameter,
    AssignmentRule,
    FunctionDefinition,
]


def import_SBML(
    filename: str, 
    name: str = None, 
    gillespy_model: "Model" = None
) -> Tuple["Model", List[str]]:
    """
    SBML to GillesPy model converter. NOTE: non-mass-action rates
    in terms of concentrations may not be converted for population
    simulation. Use caution when importing SBML.

    :param filename: Path to the SBML file for conversion.
    :param name: Name of the resulting model
    :param gillespy_model: If desired, the SBML model may be added to an existing GillesPy model

    :returns: A tuple which contains (1) the converted GillesPy2 Model and (2) a list of error strings.

    :raises ImportError: If :class:`gillespy2.sbml.SBMLImport` could not be imported.
    """

    try:
        from gillespy2.sbml.SBMLimport import convert
    except ImportError:
        raise ImportError('SBML conversion not imported successfully')

    return convert(filename, model_name=name, gillespy_model=gillespy_model)


def export_SBML(gillespy_model: "Model", filename: str = None) -> str:
    """
    GillesPy model to SBML converter

    :param gillespy_model: GillesPy model to be converted to SBML
    :param filename: Path to the SBML file for conversion. If set to None then the SBML document
        will be exported to a file with name :code:`f"{model.name}.xml"`.

    :returns: The path of exported SBML document.

    :raises ImportError: If :class:`gillespy2.sbml.SBMLExport` could not be imported.
    """
    try:
        from gillespy2.sbml.SBMLexport import export
    except ImportError:
        raise ImportError('SBML export conversion not imported successfully')

    return export(gillespy_model, path=filename)


def export_StochSS(
    gillespy_model: "Model", 
    filename: str = None, 
    return_stochss_model: bool = False
) -> Union[str, Dict]:
    """
    GillesPy model to StochSS converter

    :param gillespy_model: GillesPy model to be converted to StochSS
    :param filename: Path to the StochSS file for conversion

    :returns: If `return_stochss_model` is True then a this function will return a StochSS model dictionary. 
        If False, then the path of the StochSS model file is returned instead.

    :raises ImportError: If :class:`gillespy2.stochss.StochSSExport` could not be imported. 
    """
    try:
        from gillespy2.stochss.StochSSexport import export
    except ImportError:
        raise ImportError('StochSS export conversion not imported successfully')

    return export(gillespy_model, path=filename, return_stochss_model=return_stochss_model)


class Model(SortableObject, Jsonify):
    """
    Representation of a well mixed biochemical model. Contains reactions,
    parameters, species.

    :param name: The name of the model, or an annotation describing it.
    
    :param population: The type of model being described. A discrete stochastic model is a
        population model (True), a deterministic model is a concentration model
        (False). Automatic conversion from population to concentration models
        may be used, by setting the volume parameter.

    :param volume: The volume of the system matters when converting to from population to
        concentration form. This will also set a parameter "vol" for use in
        custom (i.e. non-mass-action) propensity functions.  
    
    :param tspan: The timepoints at which the model should be simulated. If None, a
        default timespan is added. Can also be set with :meth:`~gillespy2.core.model.Model.timespan`.
    
    :param annotation: Optional further description of model
    """

    # reserved names for model species/parameter names, volume, and operators.
    reserved_names = ['vol']
    special_characters = ['[', ']', '+', '-', '*', '/', '.', '^']

    def __init__(
        self, 
        name: str = "", 
        population: bool = True, 
        volume: float = 1.0, 
        tspan: Union[TimeSpan, Iterable[int]] = None, 
        annotation: str = "model"
    ):
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

    def __str__(self) -> str:
        divider = '\n**********\n'

        def decorate(header):
            return '\n' + divider + header + divider

        print_string = self.name
        if len(self.listOfSpecies):
            print_string += decorate('Species')
            for s in sorted(self.listOfSpecies.values()):
                print_string += '\n' + str(s)
        if len(self.listOfParameters):
            print_string += decorate('Parameters')
            for p in sorted(self.listOfParameters.values()):
                print_string += '\n' + str(p)
        if len(self.listOfReactions):
            print_string += decorate('Reactions')
            for r in sorted(self.listOfReactions.values()):
                print_string += '\n' + str(r)
        if len(self.listOfEvents):
            print_string += decorate('Events')
            for e in sorted(self.listOfEvents.values()):
                print_string += '\n' + str(e)
        if len(self.listOfAssignmentRules):
            print_string += decorate('Assignment Rules')
            for ar in sorted(self.listOfAssignmentRules.values()):
                print_string += '\n' + str(ar)
        if len(self.listOfRateRules):
            print_string += decorate('Rate Rules')
            for rr in sorted(self.listOfRateRules.values()):
                print_string += '\n' + str(rr)
        if len(self.listOfFunctionDefinitions):
            print_string += decorate('Function Definitions')
            for fd in sorted(self.listOfFunctionDefinitions.values()):
                print_string += '\n' + str(fd)

        return print_string

    def add(self, components: Union[Component, List[Component]]) -> Union[Component, List[Component]]:
        """
        Adds a component, or list of components to the model. 

        If a list is provided, Species and Parameters are added before other components.  
        Lists may contain any combination of accepted types other than lists and do not need to be in any particular order.

        The :code:`components` argument can be one or more of the following types:
            - :class:`gillespy2.Species`
            - :class:`gillespy2.Parameter`
            - :class:`gillespy2.Reaction`
            - :class:`gillespy2.Event`
            - :class:`gillespy2.RateRule`
            - :class:`gillespy2.AssignmentRule`
            - :class:`gillespy2.FunctionDefinition`
            - :class:`gillespy2.TimeSpan`

        :param components: The component or list of components to be added the the model.
        :returns: The components that were added to the model.

        :raises ModelError: :code:`components` or an element within are of incorrect type.
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

    def make_translation_table(self) -> TranslationTable:
        from collections import ChainMap

        species = self.listOfSpecies.values()
        reactions = self.listOfReactions.values()
        parameters = self.listOfParameters.values()
        assignments = self.listOfAssignmentRules.values()
        rates = self.listOfRateRules.values()
        events = self.listOfEvents.values()
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

    def serialize(self) -> str:
        """ Serializes the Model object to valid StochML. """
        self.resolve_parameters()
        doc = StochMLDocument().from_model(self)
        return doc.to_string()

    def update_namespace(self) -> None:
        """ Create a dict with flattened parameter and species objects. """
        self.namespace = OrderedDict([])
        for param in self.listOfParameters:
            self.namespace[param] = self.listOfParameters[param].value

    def sanitized_species_names(self) -> Dict[str, str]:
        """
        Generate a dictionary mapping user chosen species names to simplified formats which will be used
        later on by GillesPySolvers evaluating reaction propensity functions.

        :returns: A dictionary which maps user-defined :class:`gillespy2.Species` names to their 
            internal gillespy2 notation.
        """
        species_name_mapping = OrderedDict([])
        for i, name in enumerate(self.listOfSpecies.keys()):
            species_name_mapping[name] = 'S[{}]'.format(i)
        return species_name_mapping

    def problem_with_name(self, name: str) -> None:
        if name in Model.reserved_names:
            raise ModelError(
                'Name "{}" is unavailable. It is reserved for internal GillesPy use. Reserved Names: ({}).'.format(name,
                                                                                                                   Model.reserved_names))
        if name in self.listOfSpecies:
            raise ModelError('Name "{}" is unavailable. A species with that name exists.'.format(name))
        if name in self.listOfParameters:
            raise ModelError('Name "{}" is unavailable. A parameter with that name exists.'.format(name))
        if name in self.listOfReactions:
            raise ModelError('Name "{}" is unavailable. A reaction with that name exists.'.format(name))
        if name in self.listOfEvents:
            raise ModelError('Name "{}" is unavailable. An event with that name exists.'.format(name))
        if name in self.listOfRateRules:
            raise ModelError('Name "{}" is unavailable. A rate rule with that name exists.'.format(name))
        if name in self.listOfAssignmentRules:
            raise ModelError('Name "{}" is unavailable. An assignment rule with that name exists.'.format(name))
        if name in self.listOfFunctionDefinitions:
            raise ModelError('Name "{}" is unavailable. A function definition with that name exists.'.format(name))
        if name.isdigit():
            raise ModelError('Name "{}" is unavailable. Names must not be numeric strings.'.format(name))
        for special_character in Model.special_characters:
            if special_character in name:
                raise ModelError(
                    'Name "{}" is unavailable. Names must not contain special characters: {}.'.format(name,
                                                                                                      Model.special_characters))

    def get_species(self, s_name: str) -> Species:
        """
        Returns a species object by name.

        :param s_name: Name of the species object to be returned:
        """

        return self.listOfSpecies[s_name]

    def get_all_species(self) -> Dict[str, Species]:
        """
        :returns: A dict of all species in the model, of the form:
            {name : species object}
        """
        return self.listOfSpecies

    def add_species(self, obj: Union[Species, List[Species]]) -> Union[Species, List[Species]]:
        """
        Adds a species, or list of species to the model.

        :param obj: The species or list of species to be added to the model object
        """

        if isinstance(obj, list):
            for S in sorted(obj):
                self.add_species(S)
        else:
            try:
                self.problem_with_name(obj.name)
                self.listOfSpecies[obj.name] = obj
                self._listOfSpecies[obj.name] = 'S{}'.format(len(self._listOfSpecies))
            except Exception as e:
                raise ParameterError("Error using {} as a Species. Reason given: {}".format(obj, e))
        return obj

    def delete_species(self, name: str) -> None:
        """
        Removes a species object by name.

        :param name: Name of the species object to be removed
        """
        self.listOfSpecies.pop(name)
        if name in self._listOfSpecies:
            self._listOfSpecies.pop(name)

    def delete_all_species(self) -> None:
        """
        Removes all species from the model object.
        """
        self.listOfSpecies.clear()
        self._listOfSpecies.clear()

    def set_units(self, units: str) -> None:
        """
        Sets the units of the model to either "population" or "concentration"

        :param units: Either "population" or "concentration"

        :raises ModelError: If :code:`units.lower()` is some other value.
        """
        if units.lower() == 'concentration' or units.lower() == 'population':
            self.units = units.lower()
        else:
            raise ModelError("units must be either concentration or population (case insensitive)")

    def sanitized_parameter_names(self) -> Dict[str, str]:
        """
        Generate a dictionary mapping user chosen parameter names to simplified formats which will be used
        later on by GillesPySolvers evaluating reaction propensity functions.

        :returns: the dictionary mapping user parameter names to their internal GillesPy notation.
        """
        parameter_name_mapping = OrderedDict()
        parameter_name_mapping['vol'] = 'V'
        for i, name in enumerate(self.listOfParameters.keys()):
            if name not in parameter_name_mapping:
                parameter_name_mapping[name] = 'P{}'.format(i)
        return parameter_name_mapping

    def get_parameter(self, p_name: str) -> Parameter:
        """
        Returns a parameter object by name.

        :param p_name: Name of the parameter object to be returned
        """
        try:
            return self.listOfParameters[p_name]
        except:
            raise ModelError("No parameter named " + p_name)

    def get_all_parameters(self) -> Dict[str, Parameter]:
        """
        :returns: A dict of all parameters in the model, of the form:
            {name : parameter object}
        """
        return self.listOfParameters

    def add_parameter(self, params: Union[Parameter, List[Parameter]]) -> Union[Parameter, List[Parameter]]:
        """
        Adds a parameter, or list of parameters to the model.

        :param params:  The parameter or list of parameters to be added to the model object.
        """
        if isinstance(params, list):
            for p in sorted(params):
                self.add_parameter(p)
        else:
            if isinstance(params, Parameter) or type(params).__name__ == 'Parameter':
                self.problem_with_name(params.name)
                self.update_namespace()
                params._evaluate(self.namespace)
                self.listOfParameters[params.name] = params
                self._listOfParameters[params.name] = 'P{}'.format(len(self._listOfParameters))
            else:
                raise ParameterError("Parameter {}  must be of type {}, it is of type {}".format(params, str(type(Parameter)), str(params) ))
        return params

    def delete_parameter(self, name: str) -> None:
        """
        Removes a parameter object by name.

        :param name: Name of the parameter object to be removed
        """
        self.listOfParameters.pop(name)
        if name in self._listOfParameters:
            self._listOfParameters.pop(name)

    def set_parameter(self, p_name: str, expression: str) -> None:
        """
        Set the value of an existing parameter "pname" to "expression".

        :param p_name: Name of the parameter whose value will be set.

        :param expression: String that may be executed in C, describing the value of the
            parameter. May reference other parameters by name. (e.g. "k1*4")
        """

        p = self.listOfParameters[p_name]
        p.expression = expression
        p._evaluate()

    def resolve_parameters(self) -> None:
        """ Internal function:
        attempt to resolve all parameter expressions to scalar floats.
        This methods must be called before exporting the model.
        """
        self.update_namespace()
        for param in self.listOfParameters:
            self.listOfParameters[param]._evaluate(self.namespace)

    def delete_all_parameters(self) -> None:
        """ Deletes all parameters from model. """
        self.listOfParameters.clear()
        self._listOfParameters.clear()

    def validate_reactants_and_products(self, reactions: List[Reaction]) -> None:
        for reactant in list(reactions.reactants.keys()):
            if isinstance(reactant, str):
                if reactant not in self.listOfSpecies.keys():
                    raise ModelError(
                        'reactant: {0} for reaction {1} -- not found in model.listOfSpecies'.format(reactant,
                                                                                                    reactions.name))
                reactions.reactants[self.listOfSpecies[reactant]] = reactions.reactants[reactant]
                del reactions.reactants[reactant]
        for product in list(reactions.products.keys()):
            if isinstance(product, str):
                if product not in self.listOfSpecies.keys():
                    raise ModelError('product: {0} for reaction {1} -- not found in model.listOfSpecies'.format(product,
                                                                                                                reactions.name))
                reactions.products[self.listOfSpecies[product]] = reactions.products[product]
                del reactions.products[product]

    def add_reaction(self, reactions: Union[Reaction, List[Reaction]]) -> Union[Reaction, List[Reaction]]:
        """
        Adds a reaction, or list of reactions to the model.

        :param reactions: The reaction or list of reaction objects to be added to the model object.

        :raises ModelError: If :meth:`~gillespy2.core.Model.problem_with_name` invalidates the name of one or more 
            :class:`~gillespy2.core.reaction Reaction` objects.
        :raises ModelError: If :meth:`~gillespy2.core.Model.validate_reactants_and_products` invalidates the reactants
            and/or products on a :class:`~gillespy2.core.reaction.Reaction` object.
        :raises ModelError: If :code:`reactions` contains multiple :class:`~gillespy2.core.reaction.Reaction` objects with the same name.

        :raises ReactionError: If :meth:`gillespy2.core.reaction.Reaction.verify` invalidates the
            :class:`~gillespy2.core.reaction.Reaction` object.

        :raises ParameterError: If some other error occurs while adding the :class:`~gillespy2.core.reaction.Reaction` object.
        """

        # TODO, make sure that you cannot overwrite an existing reaction
        if isinstance(reactions, list):
            for r in sorted(reactions):
                self.add_reaction(r)
        else:
            try:
                self.problem_with_name(reactions.name)
                reactions.verify()
                self.validate_reactants_and_products(reactions)
                if reactions.name is None or reactions.name == '':
                    i = 0
                    while True:
                        if 'reaction{}'.format(i) in self.listOfReactions:
                            i += 1
                        else:
                            self.listOfReactions['reaction{}'.format(i)] = reactions
                            break
                else:
                    if reactions.name in self.listOfReactions:
                        raise ModelError("Duplicate name of reaction: {0}".format(reactions.name))
                    self.listOfReactions[reactions.name] = reactions
                # Build Sanitized reaction as well
                sanitized_reaction = Reaction(name='R{}'.format(len(self._listOfReactions)))
                sanitized_reaction.reactants = {self._listOfSpecies[species.name]: reactions.reactants[species] for
                                                species in reactions.reactants}
                sanitized_reaction.products = {self._listOfSpecies[species.name]: reactions.products[species] for
                                               species in reactions.products}
                sanitized_reaction.propensity_function = reactions.sanitized_propensity_function(self._listOfSpecies,
                                                                                                 self._listOfParameters)
                self._listOfReactions[reactions.name] = sanitized_reaction
            except Exception as e:
                raise ParameterError("Error using {} as a Reaction. Reason given: {}".format(reactions, e))
        return reactions

    def add_rate_rule(self, rate_rules: Union[RateRule, List[RateRule]]) -> Union[RateRule, List[RateRule]]:
        """
        Adds a rate rule, or list of rate rules to the model.

        :param rate_rules: The rate rule or list of rate rule objects to be added to the model object.
        """
        if isinstance(rate_rules, list):
            for rr in sorted(rate_rules):
                self.add_rate_rule(rr)
        else:
            try:
                self.problem_with_name(rate_rules.name)
                if len(self.listOfAssignmentRules) != 0:
                    for i in self.listOfAssignmentRules.values():
                        if rate_rules.variable == i.variable:
                            raise ModelError("Duplicate variable in rate_rules AND assignment_rules: {0}".
                                             format(rate_rules.variable))
                for i in self.listOfRateRules.values():
                    if rate_rules.variable == i.variable:
                        raise ModelError("Duplicate variable in rate_rules: {0}".format(rate_rules.variable))
                if rate_rules.name in self.listOfRateRules:
                    raise ModelError("Duplicate name of rate_rule: {0}".format(rate_rules.name))
                if rate_rules.formula == '':
                    raise ModelError('Invalid Rate Rule. Expression must be a non-empty string value')
                if rate_rules.variable == None:
                    raise ModelError('A GillesPy2 Rate Rule must be associated with a valid variable')
                if isinstance(rate_rules.variable, str):
                    v = rate_rules.variable
                    if v not in self.listOfSpecies and v not in self.listOfParameters:
                        raise ModelError(
                            'Invalid variable entered for Rate Rule: {}'.format(rate_rules.name))

                self.listOfRateRules[rate_rules.name] = rate_rules
                sanitized_rate_rule = RateRule(name='RR{}'.format(len(self._listOfRateRules)))
                sanitized_rate_rule.formula = rate_rules.sanitized_formula(self._listOfSpecies,
                                                                           self._listOfParameters)
                self._listOfRateRules[rate_rules.name] = sanitized_rate_rule
            except Exception as e:
                raise ParameterError("Error using {} as a Rate Rule. Reason given: {}".format(rate_rules, e))
        return rate_rules

    def add_event(self, event: Union[Event, List[Event]]) -> Union[Event, List[Event]]:
        """
        Adds an event, or list of events to the model.

        :param event: The event or list of event objects to be added to the model object.
        :type event: Event, or list of Events
        """

        if isinstance(event, list):
            for e in event:
                self.add_event(e)
        else:
            try:
                self.problem_with_name(event.name)
                if event.trigger is None or not hasattr(event.trigger, 'expression'):
                    raise ModelError(
                        'An Event must contain a valid trigger.')
                for a in event.assignments:
                    if isinstance(a.variable, str):
                        a.variable = self.get_element(a.variable)
                self.listOfEvents[event.name] = event
            except Exception as e:
                raise ParameterError("Error using {} as Event. Reason given: {}".format(event, e))
        return event

    def add_function_definition(self, function_definitions: Union[FunctionDefinition, List[FunctionDefinition]]):
        """
        Add FunctionDefinition or list of FunctionDefinitions

        :param function_definitions: The FunctionDefinition, or list of FunctionDefinitions to be added to the model
            object.
        :type function_definitions: FunctionDefinition or list of FunctionDefinitions.
        """
        if isinstance(function_definitions, list):
            for fd in function_definitions:
                self.add_function_definition(fd)
        else:
            try:
                self.problem_with_name(function_definitions.name)
                self.listOfFunctionDefinitions[function_definitions.name] = function_definitions
            except Exception as e:
                raise ParameterError(
                    "Error using {} as a Function Definition. Reason given: {}".format(function_definitions, e))

    def add_assignment_rule(self, assignment_rules: Union[AssignmentRule, List[AssignmentRule]]):
        """
        Add AssignmentRule or list of AssignmentRules to the model object.

        :param assignment_rules: The AssignmentRule or list of AssignmentRules to be added to the model object.
        :type assignment_rules: AssignmentRule or list of AssignmentRules
        """
        if isinstance(assignment_rules, list):
            for ar in assignment_rules:
                self.add_assignment_rule(ar)
        else:
            try:
                self.problem_with_name(assignment_rules.name)
                if len(self.listOfRateRules) != 0:
                    for i in self.listOfRateRules.values():
                        if assignment_rules.variable == i.variable:
                            raise ModelError("Duplicate variable in rate_rules AND assignment_rules: {0}".
                                             format(assignment_rules.variable))
                for i in self.listOfAssignmentRules.values():
                    if assignment_rules.variable == i.variable:
                        raise ModelError("Duplicate variable in assignment_rules: {0}"
                                         .format(assignment_rules.variable))
                if assignment_rules.name in self.listOfAssignmentRules:
                    raise ModelError("Duplicate name in assignment_rules: {0}".format(assignment_rules.name))
                if assignment_rules.formula == '':
                    raise ModelError('Invalid Assignment Rule. Expression must be a non-empty string value')
                if assignment_rules.variable == None:
                    raise ModelError('A GillesPy2 Rate Rule must be associated with a valid variable')

                self.listOfAssignmentRules[assignment_rules.name] = assignment_rules
            except Exception as e:
                raise ParameterError("Error using {} as a Assignment Rule. Reason given: {}".format(assignment_rules, e))

    def timespan(self, time_span: Union[TimeSpan, Iterable[int]]):
        """
        Set the time span of simulation. StochKit does not support non-uniform
        timespans. 

        :param time_span: Evenly-spaced list of times at which to sample the species populations during the simulation. 
            Best to use the form gillespy2.TimeSpan(np.linspace(<start time>, <end time>, <number of time-points, inclusive>))
        :type time_span: gillespy2.TimeSpan | iterator
        """        
        if isinstance(time_span, TimeSpan) or type(time_span).__name__ == "TimeSpan":
            self.tspan = time_span
        else:
            self.tspan = TimeSpan(time_span)

    def get_reaction(self, rname: str) -> Reaction:
        """
        :param rname: name of reaction to return
        :returns: Reaction object
        """
        return self.listOfReactions[rname]

    def get_all_reactions(self) -> Dict[str, Reaction]:
        """
        :returns: dict of all Reaction objects
        """
        return self.listOfReactions

    def delete_reaction(self, name: str) -> None:
        """
        Removes a reaction object by name.

        :param name: Name of Reaction to be removed,
        :type name: str
        """
        self.listOfReactions.pop(name)
        if name in self._listOfReactions:
            self._listOfReactions.pop(name)

    def delete_all_reactions(self) -> None:
        """
        Clears all reactions in model
        """
        self.listOfReactions.clear()
        self._listOfReactions.clear()

    def get_event(self, ename: str) -> Event:
        """
        :param ename: Name of Event to get
        :returns: Event object
        """
        return self.listOfEvents[ename]

    def get_all_events(self) -> Dict[str, Event]:
        """
        :returns: dict of all Event objects
        """
        return self.listOfEvents

    def delete_event(self, name: str) -> None:
        """
        Removes specified Event from model

        :param name: Name of Event to be removed.
        :type name: str
        """
        self.listOfEvents.pop(name)
        if name in self._listOfEvents:
            self._listOfEvents.pop(name)

    def delete_all_events(self) -> None:
        """
        Clears models events
        """
        self.listOfEvents.clear()
        self._listOfEvents.clear()

    def get_rate_rule(self, rname: str) -> RateRule:
        """
        :param rname: Name of Rate Rule to get
        :returns: RateRule object
        """
        return self.listOfRateRules[rname]

    def get_all_rate_rules(self) -> Dict[str, RateRule]:
        """
        :returns: dict of all Rate Rule objects
        """
        return self.listOfRateRules

    def delete_rate_rule(self, name: str) -> None:
        """
        Removes specified Rate Rule from model

        :param name: Name of Rate Rule to be removed.
        :type name: str
        """
        self.listOfRateRules.pop(name)
        if name in self._listOfRateRules:
            self._listOfRateRules.pop(name)

    def delete_all_rate_rules(self) -> None:
        """
        Clears all of models Rate Rules
        """
        self.listOfRateRules.clear()
        self._listOfRateRules.clear()

    def get_assignment_rule(self, aname: str) -> AssignmentRule:
        """
        :param aname: Name of Assignment Rule to get
        :returns: Assignment Rule object
        """
        return self.listOfAssignmentRules[aname]

    def get_all_assignment_rules(self) -> Dict[str, AssignmentRule]:
        """
        :returns: dict of models Assignment Rules
        """
        return self.listOfAssignmentRules

    def delete_assignment_rule(self, name: str) -> None:
        """
        Removes an assignment rule from a model

        :param name: Name of AssignmentRule object to be removed from model.
        :type name: str
        """
        self.listOfAssignmentRules.pop(name)
        if name in self._listOfAssignmentRules:
            self._listOfAssignmentRules.pop(name)

    def delete_all_assignment_rules(self) -> None:
        """
        Clears all assignment rules from model
        """
        self.listOfAssignmentRules.clear()
        self._listOfAssignmentRules.clear()

    def get_function_definition(self, fname: str) -> FunctionDefinition:
        """
        :param fname: name of Function to get
        :returns: FunctionDefinition object
        """
        return self.listOfFunctionDefinitions[fname]

    def get_all_function_definitions(self) -> Dict[str, FunctionDefinition]:
        """
        :returns: Dict of models function definitions
        """
        return self.listOfFunctionDefinitions

    def delete_function_definition(self, name: str) -> None:
        """
        Removes specified Function Definition from model

        :param name: Name of Function Definition to be removed
        :type name: str
        """
        self.listOfFunctionDefinitions.pop(name)
        if name in self._listOfFunctionDefinitions:
            self._listOfFunctionDefinitions.pop(name)

    def delete_all_function_definitions(self) -> None:
        """
        Clears all Function Definitions from a model
        """
        self.listOfFunctionDefinitions.clear()
        self._listOfFunctionDefinitions.clear()

    def get_element(self, ename: str) -> Union[Component, FunctionDefinition]:
        """
        Get element specified by name.

        :param ename: name of element to search for
        :returns: value of element, or 'element not found'
        """
        if ename in self.listOfReactions:
            return self.get_reaction(ename)
        if ename in self.listOfSpecies:
            return self.get_species(ename)
        if ename in self.listOfParameters:
            return self.get_parameter(ename)
        if ename in self.listOfEvents:
            return self.get_event(ename)
        if ename in self.listOfRateRules:
            return self.get_rate_rule(ename)
        if ename in self.listOfAssignmentRules:
            return self.get_assignment_rule(ename)
        if ename in self.listOfFunctionDefinitions:
            return self.get_function_definition(ename)
        raise ModelError(f"model.get_element(): element={ename} not found")


    def get_best_solver(self) -> GillesPySolver:
        """
        Finds best solver for the users simulation. Currently, AssignmentRules, RateRules, FunctionDefinitions,
        Events, and Species with a dynamic, or continuous population must use the TauHybridSolver.

        :param precompile: If True, and the model contains no AssignmentRules, RateRules, FunctionDefinitions, Events,
            or Species with a dynamic or continuous population, it will choose SSACSolver
        :type precompile: bool

        :returns: gillespy2.gillespySolver
        """
        from gillespy2.solvers.numpy import can_use_numpy
        hybrid_check = False
        chybrid_check = True
        if len(self.get_all_rate_rules())  or len(self.get_all_events()):
            hybrid_check = True
        if len(self.get_all_assignment_rules()) or len(self.get_all_function_definitions()):
            hybrid_check = True
            chybrid_check = False

        if len(self.get_all_species()) and hybrid_check == False:
            for i in self.get_all_species():
                tempMode = self.get_species(i).mode
                if tempMode == 'dynamic' or tempMode == 'continuous':
                    hybrid_check = True
                    break

        from gillespy2.solvers.cpp.build.build_engine import BuildEngine
        can_use_cpp = not len(BuildEngine.get_missing_dependencies())

        if not can_use_cpp and not can_use_numpy:
            raise ModelError('Dependency Error, cannot run model.')

        if can_use_cpp and hybrid_check and chybrid_check:
            from gillespy2 import TauHybridCSolver
            return TauHybridCSolver
        elif can_use_numpy and hybrid_check:
            from gillespy2 import TauHybridSolver
            return TauHybridSolver
        
        if can_use_cpp is False and can_use_numpy and not hybrid_check:
            from gillespy2 import NumPySSASolver
            return NumPySSASolver

        else:
            from gillespy2 import SSACSolver
            return SSACSolver

    def get_best_solver_algo(self, algorithm: str) -> GillesPySolver:
        """
        If user has specified a particular algorithm, we return either the Python or C++ version of that algorithm
        """
        from gillespy2.solvers.numpy import can_use_numpy
        from gillespy2.solvers.cpp.build.build_engine import BuildEngine
        can_use_cpp = not len(BuildEngine.get_missing_dependencies())
        chybrid_check = True
        if len(self.get_all_assignment_rules()) or len(self.get_all_function_definitions()):
            chybrid_check = False

        if not can_use_cpp and can_use_numpy:
            raise ModelError("Please install C++ or Numpy to use GillesPy2 solvers.")

        if algorithm == 'Tau-Leaping':
            if can_use_cpp:
                from gillespy2 import TauLeapingCSolver
                return TauLeapingCSolver
            else:
                from gillespy2 import TauLeapingSolver
                return TauLeapingSolver

        elif algorithm == 'SSA':
            if can_use_cpp:
                from gillespy2 import SSACSolver
                return SSACSolver
            else:
                from gillespy2 import NumPySSASolver
                return NumPySSASolver

        elif algorithm == 'ODE':
            if can_use_cpp:
                from gillespy2 import ODECSolver
                return ODECSolver
            else:
                from gillespy2 import ODESolver
                return ODESolver

        elif algorithm == 'Tau-Hybrid':
            if can_use_cpp and chybrid_check:
                from gillespy2 import TauHybridCSolver
                return TauHybridCSolver
            else:
                from gillespy2 import TauHybridSolver
                return TauHybridSolver

        elif algorithm == 'CLE':
            from gillespy2 import CLESolver
            return CLESolver
            
        else:
            raise ModelError("Invalid value for the argument 'algorithm' entered. "
                             "Please enter 'SSA', 'ODE', 'CLE', 'Tau-leaping', or 'Tau-Hybrid'.")

    def get_model_features(self) -> AbstractSet[Union[Event, RateRule, AssignmentRule, FunctionDefinition]]:
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

    def run(
        self, 
        solver: Optional[GillesPySolver] = None, 
        timeout: int = 0, 
        t: int = None, 
        increment: float = None, 
        show_labels: bool = True, 
        algorithm: str = None,
        **solver_args
    ) -> Results:
        """
        Run a simulation with the current :class:`~gillespy2.core.model.Model` instance.

        A solver can be specified to tweak simulation results and improve simulation runtime. This package offers 
        multiple different algorithms implemented in either C++ or Python. We recommend the use of 
        :meth:`~gillespy2.core.model.Model.get_best_solver` to determine the best solver for your use-case.
        
        Simulations can be paused by issuing a keyboard interrupt during runtime. This is achieved by pressing 
        Control + C or the *stop* button within a Jupyter Notebook session. The :class:`~gillespy2.core.results.Results`
        object returned after the interrupt can be used to resume the simulation at a later time.

        To resume a simulation, pass the :class:`~gillespy2.core.results.Results` object to this method,
        ensuring that :code:`t` is set to the time the resumed simulation will end.

        .. code-block:: python

            model = MichaelisMenten()
            results = model.run()

            # Control + C 

            results = model.run(resume=results, t=x)

        .. warning::
            Pause / resume functionality is only supported for single trajectory simulations. :code:`t` must be set or 
            unexpected behavior may occur.

            .. code-block:: text

                <>=======()
               (/\___   /|\\\\          ()==========<>_
                     \_/ | \\\\        //|\   ______/ \)
                       \_|  \\\\      // | \_/
                         \|\/|\_   //  /\/
                          (--)\ \_//  /
                         //_/\_\/ /  |
                        @@/  |=\  \  |
                             \_=\_ \ |
                               \==\ \|\_ 
                            __(\===\(  )\\
                           (((~) __(_/   |
                                (((~) \  /
                                ______/ /
                                '------'

        .. note::
            :code:`show_labels = False` has been deprecated and may not be supported in future versions of this software.

        :param solver: The solver by which to simulate the current :class:`~gillespy2.core.model.Model`.
            An already initialized solver can also be passed via this keyword argument. A solver will be automatically
            chosen if this argument is :code:`None`. 

            See :code:`algorithm` to specify an algorithm *type*.

        :param timeout: The timeout, in seconds, to restrict simulation runtime.

        :param t: The end time of the simulation.

        :param algorithm: The algorithm to compute the simulation with. A specific solver will be automatically chosen
            based on the value of this argument with :meth:`~gillespy2.core.model.Model.get_best_solver_algo`.

            Possible values: :code:`"ODE"`, :code:`"Tau-Leaping"`, or :code:`"SSA"`.

        :param solver_args: Additional solver-specific arguments to be passed to :code:`solver.run()`.

        :returns: A :class:`~gillespy2.core.results.Results` object which contains one or more computed
            :class:`~gillespy2.core.trajectory.Trajectory` objects. Helper functions on the results object can
            then be used to graph, manipulate, or export the data.

            .. note::
                If :code:`show_labels = False` this function will return a two-dimensional :class:`numpy.ndarray`
                that contains species population data.

        """

        if not show_labels:
            from gillespy2.core import log
            log.warning('show_labels = False is deprecated. Future releases of GillesPy2 may not support this feature.')

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
        except Exception as e:
            raise SimulationError(
                "argument 'solver={}' to run() failed.  Reason Given: {}".format(solver, e)
            ) from e


class StochMLDocument():
    """ Serializiation and deserialization of a Model to/from
        the native StochKit2 XML format. """

    def __init__(self):
        # The root element
        self.document = eTree.Element("Model")
        self.annotation = None

    @classmethod
    def from_model(cls, model: Model) -> "StochMLDocument":
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
        md = cls()

        d = eTree.Element('Description')

        #
        if model.units.lower() == "concentration":
            d.set('units', model.units.lower())

        d.text = model.annotation
        md.document.append(d)

        # Number of Reactions
        nr = eTree.Element('NumberOfReactions')
        nr.text = str(len(model.listOfReactions))
        md.document.append(nr)

        # Number of Species
        ns = eTree.Element('NumberOfSpecies')
        ns.text = str(len(model.listOfSpecies))
        md.document.append(ns)

        # Species
        spec = eTree.Element('SpeciesList')
        for sname in model.listOfSpecies:
            spec.append(md.__species_to_element(model.listOfSpecies[sname]))
        md.document.append(spec)

        # Parameters
        params = eTree.Element('ParametersList')
        for pname in model.listOfParameters:
            params.append(md.__parameter_to_element(
                model.listOfParameters[pname]))

        params.append(md.__parameter_to_element(Parameter(name='vol', expression=model.volume)))

        md.document.append(params)

        # Reactions
        reacs = eTree.Element('ReactionsList')
        for rname in model.listOfReactions:
            reacs.append(md.__reaction_to_element(model.listOfReactions[rname], model.volume))
        md.document.append(reacs)

        return md

    @classmethod
    def from_file(cls, filepath: str) -> "StochMLDocument":
        """ Intializes the document from an exisiting native StochKit XML
        file read from disk. """
        tree = eTree.parse(filepath)
        root = tree.getroot()
        md = cls()
        md.document = root
        return md

    @classmethod
    def from_string(cls, string: str) -> "StochMLDocument":
        """ Intializes the document from an exisiting native StochKit XML
        file read from disk. """
        root = eTree.fromString(string)

        md = cls()
        md.document = root
        return md

    def to_model(self, name: str) -> Model:
        """ Instantiates a Model object from a StochMLDocument. """

        # Empty model
        model = Model(name=name)
        root = self.document

        # Try to set name from document
        if model.name == "":
            name = root.find('Name')
            if name.text is None:
                raise NameError("The Name cannot be none")
            else:
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
        for px in root.iter('Parameter'):
            name = px.find('Id').text
            expr = px.find('Expression').text
            if name.lower() == 'vol' or name.lower() == 'volume':
                model.volume = float(expr)
            else:
                p = Parameter(name, expression=expr)
                # Try to evaluate the expression in the empty namespace
                # (if the expr is a scalar value)
                p._evaluate()
                model.add_parameter(p)

        # Create species
        for spec in root.iter('Species'):
            name = spec.find('Id').text
            val = spec.find('InitialPopulation').text
            if '.' in val:
                val = float(val)
            else:
                val = int(val)
            s = Species(name, initial_value=val)
            model.add_species([s])

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
            except:
                raise InvalidStochMLError("Reaction has no name.")

            reaction = Reaction(name=name, reactants={}, products={})

            # Type may be 'mass-action','customized'
            try:
                type = reac.find('Type').text
            except:
                raise InvalidStochMLError("No reaction type specified.")

            reactants = reac.find('Reactants')
            try:
                for ss in reactants.iter('SpeciesReference'):
                    specname = ss.get('id')
                    # The stochiometry should be an integer value, but some
                    # exising StoxhKit models have them as floats. This is
                    # why we need the slightly odd conversion below.
                    stoch = int(float(ss.get('stoichiometry')))
                    # Select a reference to species with name specname
                    sref = model.listOfSpecies[specname]
                    try:
                        # The sref list should only contain one element if
                        # the XML file is valid.
                        reaction.reactants[sref] = stoch
                    except Exception as e:
                        StochMLImportError(e)
            except:
                # Yes, this is correct. 'reactants' can be None
                pass

            products = reac.find('Products')
            try:
                for ss in products.iter('SpeciesReference'):
                    specname = ss.get('id')
                    stoch = int(float(ss.get('stoichiometry')))
                    sref = model.listOfSpecies[specname]
                    try:
                        # The sref list should only contain one element if
                        # the XML file is valid.
                        reaction.products[sref] = stoch
                    except Exception as e:
                        raise StochMLImportError(e)
            except:
                # Yes, this is correct. 'products' can be None
                pass

            if type == 'mass-action':
                reaction.massaction = True
                reaction.type = 'mass-action'
                # If it is mass-action, a parameter reference is needed.
                # This has to be a reference to a species instance. We
                # explicitly disallow a scalar value to be passed as the
                # parameter.
                try:
                    ratename = reac.find('Rate').text
                    try:
                        reaction.marate = model.listOfParameters[ratename]
                    except KeyError as k:
                        # No paramter name is given. This is a valid use case
                        # in StochKit. We generate a name for the paramter,
                        # and create a new parameter instance. The parameter's
                        # value should now be found in 'ratename'.
                        generated_rate_name = "Reaction_" + name + \
                                              "_rate_constant"
                        p = Parameter(name=generated_rate_name,
                                      expression=ratename)
                        # Try to evaluate the parameter to set its value
                        p._evaluate()
                        model.add_parameter(p)
                        reaction.marate = model.listOfParameters[
                            generated_rate_name]

                    reaction.create_mass_action()
                except Exception as e:
                    raise
            elif type == 'customized':
                try:
                    propfunc = reac.find('PropensityFunction').text
                except Exception as e:
                    raise InvalidStochMLError(
                        "Found a customized propensity function, but no expression was given. {}".format(e))
                reaction.propensity_function = propfunc
                reaction.ode_propensity_function = propfunc
            else:
                raise InvalidStochMLError(
                    "Unsupported or no reaction type given for reaction" + name)

            model.add_reaction(reaction)

        return model

    def to_string(self) -> str:
        """ Returns  the document as a string. """
        try:
            doc = eTree.tostring(self.document, pretty_print=True)
            return doc.decode("utf-8")
        except:
            # Hack to print pretty xml without pretty-print
            # (requires the lxml module).
            doc = eTree.tostring(self.document)
            xmldoc = xml.dom.minidom.parseString(doc)
            uglyXml = xmldoc.toprettyxml(indent='  ')
            text_re = re.compile(">\n\s+([^<>\s].*?)\n\s+</", re.DOTALL)
            prettyXml = text_re.sub(">\g<1></", uglyXml)
            return prettyXml

    def __species_to_element(self, S: Species):
        e = eTree.Element('Species')
        idElement = eTree.Element('Id')
        idElement.text = S.name
        e.append(idElement)

        if hasattr(S, 'description'):
            descriptionElement = eTree.Element('Description')
            descriptionElement.text = S.description
            e.append(descriptionElement)

        initialPopulationElement = eTree.Element('InitialPopulation')
        initialPopulationElement.text = str(S.initial_value)
        e.append(initialPopulationElement)

        return e

    def __parameter_to_element(self, P: Parameter):
        e = eTree.Element('Parameter')
        idElement = eTree.Element('Id')
        idElement.text = P.name
        e.append(idElement)
        expressionElement = eTree.Element('Expression')
        expressionElement.text = str(P.value)
        e.append(expressionElement)
        return e

    def __reaction_to_element(self, R: Reaction, model_volume: float):
        e = eTree.Element('Reaction')

        idElement = eTree.Element('Id')
        idElement.text = R.name
        e.append(idElement)

        descriptionElement = eTree.Element('Description')
        descriptionElement.text = self.annotation
        e.append(descriptionElement)

        # StochKit2 wants a rate for mass-action propensites
        if R.massaction and model_volume == 1.0:
            rateElement = eTree.Element('Rate')
            # A mass-action reactions should only have one parameter
            rateElement.text = R.marate.name
            typeElement = eTree.Element('Type')
            typeElement.text = 'mass-action'
            e.append(typeElement)
            e.append(rateElement)

        else:
            typeElement = eTree.Element('Type')
            typeElement.text = 'customized'
            e.append(typeElement)
            functionElement = eTree.Element('PropensityFunction')
            functionElement.text = R.propensity_function
            e.append(functionElement)

        reactants = eTree.Element('Reactants')

        for reactant, stoichiometry in R.reactants.items():
            srElement = eTree.Element('SpeciesReference')
            srElement.set('id', str(reactant.name))
            srElement.set('stoichiometry', str(stoichiometry))
            reactants.append(srElement)

        e.append(reactants)

        products = eTree.Element('Products')
        for product, stoichiometry in R.products.items():
            srElement = eTree.Element('SpeciesReference')
            srElement.set('id', str(product.name))
            srElement.set('stoichiometry', str(stoichiometry))
            products.append(srElement)
        e.append(products)

        return e
