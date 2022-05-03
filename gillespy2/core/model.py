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
    SpeciesError,
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
    SBML to gillespy2 :class:`Model` converter.

    .. note::
        Mass-action rates *not* in terms of concentration may not be converted correctly to a population-based
        simulation. Use caution when using this function to convert SBML to :class:`Model`.

    :param filename: The path of the SBML file to convert.
    :param name: The name to give the resulting :class:`Model`.
    :param gillespy2_model: An already existing :class:`Model` object that the SBML model will be added to.

    :returns: A tuple which contains (1) the converted GillesPy2 Model and (2) a list of error strings.
    :rtype: Tuple[Model, List[str]]

    :raises ImportError: If :class:`gillespy2.sbml.SBMLImport` could not be imported.
    """

    try:
        from gillespy2.sbml.SBMLimport import convert
    except ImportError:
        raise ImportError('SBML conversion not imported successfully')

    return convert(filename, model_name=name, gillespy_model=gillespy_model)


def export_SBML(gillespy_model: "Model", filename: str = None) -> str:
    """
    gillespy2 :class:`Model` to SBML converter.

    :param gillespy2_model: The gillespy2 :class:`Model` to convert to SBML.
    :param filename: The path of the new SBML document. If `None` then the SBML document will be written to a file with
        the name :code:`f"{model.name}".xml`.

    :returns: The path on disk of the exported SBML document.
    :rtype: str

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
    gillespy2 :class:`Model` to StochSS model converter.

    :param gillespy2_model: The gillespy2 :class:`Model` to be converted to a StochSS-compatible model.
    :param filename: The path on disk of the exported StochSS model.

    :returns: If :code:`return_stochss_model = True` then this function will return a StochSS model dictionary.
        If :code:`False` then the path of the StochSS model file is returned instead.
    :rtype: Union[str, Dict]

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

    :param name: The name of the model or an annotation which describes it.

    :param population: The type of model being described. A discrete stochastic model is a population model
        (:code:`population = True`), a deterministic model is a concentration model (:code:`population = False`). A
        population model can be automatically converted to a concentration model by setting the :code:`volume` parameter.

    :param volume: The volume of the system. This parameter is needed when converting a model from population to
        concentration. This also sets the parameter :code:`"vol"` for use in custom (i.e. non-mass-action) propensity
        functions.

    :param tspan: The timepoint at which the model should be simulated. If :code:`None` then a default timespan is
        selected. This model property can also be set with :meth:`timespan`.

    :param annotation: An optional description of the model.

    :raises ModelError: If :code:`population = False` and :code:`volume != 1.0`. This occurs because a concentration
        model implicitly accounts for volume. It is important to note that :code:`population = True` sets
        :code:`self.units = "population"`, and :code:`population = False` sets :code:`self.units = "concentation"`. The
        latter case, in combination with :code:`volume != 1.0`, is what causes this error to occur.

        .. note::
            Concentration models (e.g :code:`self.units = "concentration"`) may only be simulated deterministically.
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

    def problem_with_name(self, name: str) -> None:
        """
        Validation function which ensures that :code:`name`:

        1. Does not collide with reserved gillespy2 keywords.
        2. Is not already in use by another object in this :class:`Model` instance.
        3. Does not contain any reserved characters.

        This function will raise an exception if :code:`name` is invalidated by any of the previous cases.

        A list of reserved names can be accessed via the :const:`Model.reserved_names` property (e.g
        :code:`print(Model.reserved_names)`).

        A list of reserved characters can be accessed via the :const:`Model.special_characters` property (e.g
        :code:`print(Model.special_characters)`).

        :param name: The name to validate against the reserved keywords list and other already set identifiers.

        :raises ModelError: If :code:`name` is reserved for internal gillespy2 use.

        :raises ModelError: If :code:`name` is already in use by any one of the following types:

            - :class:`gillespy2.core.species.Species`
            - :class:`gillespy2.core.parameter.Parameter`
            - :class:`gillespy2.core.reaction.Reaction`
            - :class:`gillespy2.core.events.Event`
            - :class:`gillespy2.core.raterule.RateRule`
            - :class:`gillespy2.core.assignmentrule.AssignmentRule`
            - :class:`gillespy2.core.functiondefinition.FunctionDefinition`

        :raises ModelError: If :code:`name` is a numeric string -- e.g :code:`name.isdigit() == True`.

        :raises ModelError: If :code:`name` contains one or more reserved special characters.
        """
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
    def update_namespace(self) -> None:
        """
        Creates the internal :code:`namespace` property by flattening :code:`self.listOfParameters` into a dictionary
        of :code:`(parameter, parameter.value)` pairs.

        This function mutates the interior state of this :class:`Model` object.
        """
        self.namespace = OrderedDict([])
        for param in self.listOfParameters:
            self.namespace[param] = self.listOfParameters[param].value

    def add(self, components: Union[Component, List[Component]]) -> Union[Component, List[Component]]:
        """
        Adds one or more components to the :class:`Model`.

        If a list is provided then :class:`~gillespy2.core.species.Species` and
        :class:`~gillespy2.core.parameter.Parameter` objects are added before other components.

        The :code:`components` argument can be one or more of the following types:
            - :class:`gillespy2.Species`
            - :class:`gillespy2.Parameter`
            - :class:`gillespy2.Reaction`
            - :class:`gillespy2.Event`
            - :class:`gillespy2.RateRule`
            - :class:`gillespy2.AssignmentRule`
            - :class:`gillespy2.FunctionDefinition`
            - :class:`gillespy2.TimeSpan`

        :param components: The component(s) to add to the :class:`Model`.

        :returns: The components that were added to the model.
        :rtype: One or more objects with type defined in the docstring above.

        :raises ModelError: If :code:`components` contains an element of incorrect type.
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

    def add_species(self, species: Union[Species, List[Species]]) -> Union[Species, List[Species]]:
        """
        Adds one or more :class:`~gillespy2.core.species.Species` objects to this :class:`Model` instance.

        :param species: The species or list of species objects to be added to this model.

        :returns: One or more species objects that were successfully added to this :class:`Model`.
        :rtype: Union[Species, List[Species]]

        :raises ModelError: If :meth:`Model.problem_with_name` invalidates the name of one or more species objects.
        :raises ModelError: If :code:`species` contains an element of incorrect type or if :meth:`Species.validate` fails.
        :raises ModelError: If some other error occurs while adding the species object(s).
        """
        if isinstance(species, list):
            for spec in sorted(species):
                self.add_species(spec)
        elif isinstance(species, Species) or type(species).__name__ == "Species":
            try:
                species.validate()
                self.problem_with_name(species.name)
                self.listOfSpecies[species.name] = species
                self._listOfSpecies[species.name] = f'S{len(self._listOfSpecies)}'
            except SpeciesError as err:
                errmsg = f"Could not add species: {species.name}, Reason given: {err}"
                raise ModelError(errmsg) from err
        else:
            errmsg = f"species must be of type Species or list of Species not {type(species)}"
            raise ModelError(errmsg)
        return species

    def delete_species(self, name: str) -> None:
        """
        Removes a :class:`~gillespy2.core.species.Species` object from this :class:`Model` instance by name.

        :param name: The name of the species object to remove.
        :raises KeyError: If a species with name :code:`name` does not exist.
        """
        self.listOfSpecies.pop(name)
        if name in self._listOfSpecies:
            self._listOfSpecies.pop(name)

    def delete_all_species(self) -> None:
        """
        Removes all :class:`~gillespy2.core.species.Species` objects attached to this :class:`Model` instance.
        """
        self.listOfSpecies.clear()
        self._listOfSpecies.clear()

    def get_species(self, name: str) -> Species:
        """
        Returns a :class:`~gillespy2.core.species.Species` object by name.

        :param name: The name of the species object to return.

        :returns: The species object.
        :rtype: Species

        :raises KeyError: If a species with name :code:`name` does not exist.
        """
        if name not in self.listOfSpecies:
            raise ModelError(f"Species {name} could not be found in the model.")
        return self.listOfSpecies[name]

    def get_all_species(self) -> Dict[str, Species]:
        """
        Returns all :class:`~gillespy2.core.species.Species` objects attached to this :class:`Model` instance.

        :returns: A dictionary of species of the form :code:`{ name: Species }`
        :rtype: Dict[str, Species]
        """
        return self.listOfSpecies

    def sanitized_species_names(self) -> Dict[str, str]:
        """
        Generates a dictionary which maps user-defined :class:`~gillespy2.core.species.Species` names to the simplified
        formats used by solvers when evaluating :class:`~gillespy2.core.reaction.Reaction` propensity functions.

        :returns: A dictionary which maps user-defined :class:`~gillespy2.core.species.Species` names to their internal
            gillespy2 notation.

        :rtype: Dict[str, str]
        """
        species_name_mapping = OrderedDict([])
        for i, name in enumerate(self.listOfSpecies.keys()):
            species_name_mapping[name] = 'S[{}]'.format(i)
        return species_name_mapping

    def add_parameter(self, parameters: Union[Parameter, List[Parameter]]) -> Union[Parameter, List[Parameter]]:
        """
        Adds one or more :class:`~gillespy2.core.parameter.Parameter` objects to this :class:`Model` instance.

        :param parameters: The parameter or list of parameter objects to be added to the model object.

        :raises ModelError: If :meth:`Model.problem_with_name` invalidates the name of one or more
            :class:`~gillespy2.core.parameter.Parameter` objects.

        :raises ParameterError: If one or more elements within :code:`parameters` are not of type
            :class:`~gillespy2.core.parameter.Parameter`.
        """
        if isinstance(parameters, list):
            for param in sorted(parameters):
                self.add_parameter(param)
        elif isinstance(parameters, Parameter) or type(parameters).__name__ == 'Parameter':
            self.problem_with_name(parameters.name)
            self.resolve_parameter(parameters)
            self.listOfParameters[parameters.name] = parameters
            self._listOfParameters[parameters.name] = f'P{len(self._listOfParameters)}'
        else:
            errmsg = f"parameters must be of type Parameter or list of Parameter not {type(parameters)}"
            raise ModelError(errmsg)
        return parameters

    def delete_parameter(self, name: str) -> None:
        """
        Removes a :class:`~gillespy2.core.parameter.Parameter` object from this :class:`Model` by name.

        :param name: The name of the parameter object to remove.

        :raises KeyError: If a parameter with name :code:`name` does not exist.
        """
        self.listOfParameters.pop(name)
        if name in self._listOfParameters:
            self._listOfParameters.pop(name)

    def delete_all_parameters(self) -> None:
        """
        Removes all :class:`~gillespy2.core.parameter.Parameter` objects attached to this :class:`Model` instance.
        """
        self.listOfParameters.clear()
        self._listOfParameters.clear()

    def get_parameter(self, name: str) -> Parameter:
        """
        Returns a :class:`~gillespy2.core.parameter.Parameter` object by name.

        :param name: The name of the parameter object to return.

        :returns: The parameter object.
        :rtype: Parameter

        :raises ModelError: If a parameter with name :code:`name` does not exist.
        """
        if name not in self.listOfParameters:
            raise ModelError(f"Parameter {name} could not be found in the model.")
        return self.listOfParameters[name]

    def get_all_parameters(self) -> Dict[str, Parameter]:
        """
        Returns all :class:`~gillespy2.core.parameter.Parameter` objects attached to this :class:`Model` instance.

        :returns: A dictionary of species of the form :code:`{ name: Parameter }`
        :rtype: Dict[str, Species]
        """
        return self.listOfParameters

    def resolve_parameter(self, parameter: Parameter) -> None:
        """
        This function is for internal use.

        Attempt to resolve the given :class:`~gillespy2.core.parameter.Parameter` expressions to scalar floats.
        This method must be called **prior** to exporting this :class:`Model` instance.

        :param parameter: The target parameter to resolve.

        :raises ModelError: If the :code:`parameter` can't be resolved.
        """
        try:
            parameter.validate()
            self.update_namespace()
            parameter._evaluate(self.namespace)
        except ParameterError as err:
            raise ModelError(f"Could not add/resolve parameter: {parameter.name}, Reason given: {err}") from err

    def resolve_all_parameters(self) -> None:
        """
        This function is for internal use.

        Attempt to resolve all :class:`~gillespy2.core.parameter.Parameter` expressions to scalar floats. This method
        must be called **prior** to exporting this :class:`Model` instance.
        """
        for _, parameter in self.listOfParameters.items():
            self.resolve_parameter(parameter)

    def sanitized_parameter_names(self) -> Dict[str, str]:
        """
        Generates a dictionary which maps user-defined :class:`~gillespy2.core.parameter.Parameter` names to the
        simplified formats used by solvers when evaluating :class:`~gillespy2.core.reaction.Reaction` propensity
        functions.

        :returns: A dictionary which maps user-defined :class:`~gillespy2.core.parameter.Parameter` names to their
            internal gillespy2 notation.

        :rtype: Dict[str, str]
        """
        parameter_name_mapping = OrderedDict()
        parameter_name_mapping['vol'] = 'V'
        for i, name in enumerate(self.listOfParameters.keys()):
            if name not in parameter_name_mapping:
                parameter_name_mapping[name] = 'P{}'.format(i)
        return parameter_name_mapping

    def set_parameter(self, p_name: str, expression: str) -> None:
        """
        Sets the expression of a :class:`~gillespy2.core.parameter.Parameter` on this :class:`Model` instance (Deprecated).

        :param p_name: The name of the parameter whose :code:`expression` property will be set.

        :param expression: The expression to be applied, written as a string evaluatable by a solver during a
            simulation. reference other parameters by name -- e.g :code:`"k1 * 4`.
        """
        from gillespy2.core import log
        log.warning(
            """
            Model.set_parameter has been deprecated.  Future releases of GillesPy2 may
            not support this feature.  Parameter.expression should only be set in the constructor.
            """
        )

        parameter = self.listOfParameters[p_name]
        parameter.expression = expression
        self.resolve_parameter(parameter)

    def add_reaction(self, reactions: Union[Reaction, List[Reaction]]) -> Union[Reaction, List[Reaction]]:
        """
        Adds one or more :class:`~gillespy2.core.reaction.Reaction` objects to this :class:`Model` instance.

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
        if isinstance(reactions, list):
            for reaction in sorted(reactions):
                self.add_reaction(reaction)
        elif isinstance(reactions, Reaction) or type(reactions).__name__ == "Reaction":
            self.problem_with_name(reactions.name)
            self.resolve_reaction(reactions)
            self.listOfReactions[reactions.name] = reactions
            # Build Sanitized reaction as well
            sanitized_reaction = reactions._create_sanitized_reaction(
                len(self.listOfReactions), self._listOfSpecies, self._listOfParameters
            )
            self._listOfReactions[reactions.name] = sanitized_reaction
        else:
            errmsg = f"reactions must be of type Reaction or list of Reaction not {type(reactions)}"
            raise ModelError(errmsg)
        return reactions

    def delete_reaction(self, name: str) -> None:
        """
        Removes a :class:`~gillespy2.core.reaction.Reaction` object from this :class:`Model` by name.

        :param name: The name of the reaction to remove.
        :raises KeyError: If a reaction with name :code:`name` does not exist.
        """
        self.listOfReactions.pop(name)
        if name in self._listOfReactions:
            self._listOfReactions.pop(name)

    def delete_all_reactions(self) -> None:
        """
        Removes all :class:`~gillespy2.core.reaction.Reaction` objects attached to this :class:`Model` instance.
        """
        self.listOfReactions.clear()
        self._listOfReactions.clear()

    def get_reaction(self, name: str) -> Reaction:
        """
        Returns a :class:`~gillespy2.core.reaction.Reaction` on this :class:`Model` by name.

        :param name: The name of the reaction to return.

        :returns: The reaction with name :class:`name`.
        :rtype: Reaction

        :raises KeyError: If a :class:`~gillespy2.core.reaction.Reaction` with name :code:`name` does not exist on
            this :class:`Model`.
        """
        if name not in self.listOfReactions:
            raise ModelError(f"Reaction {name} could not be found in the model.")
        return self.listOfReactions[name]

    def get_all_reactions(self) -> Dict[str, Reaction]:
        """
        Returns all :class:`~gillespy2.core.reaction.Reaction` objects attached to this :class:`Model` instance.

        :returns: A dictionary of reactions of the form :code:`{ name: Reaction }`.
        :rtype: Dict[str, Reaction]
        """
        return self.listOfReactions

    def resolve_reaction(self, reaction: Reaction) -> None:
        """
        This function is for internal use.

        Validates the given :class:`~gillespy2.core.reaction.Reaction` object against this :class:`Model` instance's
        internal state.

        More specifically, this function ensures that reactants and products on a
        :class:`~gillespy2.core.reaction.Reaction` object reference :class:`~gillesp2.core.species.Species` that exist
        within :attr:`Model.listOfSpecies`.

        :param reaction: The reaction to validate.

        :raises ModelError: If a reactant in :code:`reaction` references a species in :attr:`Model.listOfSpecies` that
            does not exist.

        :raises ModelError: If a product in :code:`reaction` references a species in :attr:`Model.listOfSpecies` that
            does not exist.
        """
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

    def resolve_all_reactions(self) -> None:
        """
        This function is for internal use.

        Ensure that the rate and all reactants and products are present in the model
        for all reactions.  This methods must be called before exporting the model.
        """
        for _, reaction in self.listOfReactions.items():
            self.resolve_reaction(reaction)

    def validate_reactants_and_products(self, reactions: List[Reaction]) -> None:
        """
        This function is for internal use (Deprecated).

        Validates the given :class:`~gillespy2.core.reaction.Reaction` object against this :class:`Model` instance's
        internal state.

        More specifically, this function ensures that reactants and products on a
        :class:`~gillespy2.core.reaction.Reaction` object reference :class:`~gillesp2.core.species.Species` that exist
        within :attr:`Model.listOfSpecies`.

        :param reactions: The reactions to validate.

        :raises ModelError: If a reactant in :code:`reactions` references a species in :attr:`Model.listOfSpecies` that
            does not exist.

        :raises ModelError: If a product in :code:`reactions` references a species in :attr:`Model.listOfSpecies` that
            does not exist.
        """
        from gillespy2.core import log
        log.warning(
            """
            Model.validate_reactants_and_products has been deprecated. Future releases of
            GillesPy2 may not support this feature.  Use Model.resolve_reaction instead.
            """
        )

        self.resolve_reaction(reactions)

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
        """
        Converts the current :class:`Model` object to StochML.

        :returns: The contents of the newly created StochML document.
        :rtype: str
        """
        self.resolve_all_parameters()
        doc = StochMLDocument().from_model(self)
        return doc.to_string()

    def set_units(self, units: str) -> None:
        """
        Sets the units of this :class:`Model` instance to either :code:`"population"` or :code:`"concentration"`.

        :param units: Either :code:`"population"` or :code:`"concentration"`.

        :raises ModelError: If :code:`units.lower()` is some other value.
        """
        if units.lower() == 'concentration' or units.lower() == 'population':
            self.units = units.lower()
        else:
            raise ModelError("units must be either concentration or population (case insensitive)")

    def add_rate_rule(self, rate_rules: Union[RateRule, List[RateRule]]) -> Union[RateRule, List[RateRule]]:
        """
        Adds one or more :class:`~gillespy2.core.raterule.RateRule` objects to this :class:`Model` instance.

        :param rate_rules: The rate rule or list of rate rules to be added to this model.

        :returns: One or more rate rule objects that were successfully added to this :class:`Model`.
        :rtype: Union[RateRule, List[RateRule]]

        :raises ModelError: If :func:`Model.problem_with_name` invalidates the name of one or more
            :class:`~gillespy2.core.raterule.RateRule` objects.

        :raises ModelError: If the :code:`variable` property of a :class:`~gillespy2.core.raterule.RateRule` instance
            is already set within a :class:`~gillespy2.core.assignmentrule.AssignmentRule` or
            :class:`~gillespy2.core.raterule.RateRule` attached to this :class:`Model`.

        :raises ModelError: If one or more rate rules in :code:`rate_rules` has the same :code:`name` as a
            :class:`~gillespy2.core.raterule.RateRule` on this :class:`Model` instance.

        :raises ModelError: If one or more rate rules in :code:`rate_rules` has an invalid :attr:`RateRule.formula`
            value. Specifically, if :code:`RateRule.formula == ""`.

        :raises ModelError: If one or more rate rules in :code:`rate_rules` has a :attr:`RateRule.variable` value of
            :code:`None`.

        :raises ModelError: If one or more rate rules in :code:`rate_rules` has a :attr:`RateRule.variable` which
            references a :class:`~gillespy2.core.species.Species` or :class:`~gillespy2.core.parameter.Parameter`
            that does not exist. This function specifically checks that :attr:`RateRule.variable` exists in either
            :attr:`Model.listOfSpecies` and :attr:`Model.listOfParameters` and, if not, raises an error.

        :raises ParameterError: If some other error occurs while adding the :class:`~gillespy2.core.raterule.RateRule`
            object(s).
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
        Adds one or more :class:`~gillespy2.core.events.Event` to this :class:`Model` instance.

        :param event: The event or list of events to be added to this model.

        :returns: One or more event objects that were successfully added this :class:`Model`.
        :rtype: Union[Event, List[Event]]

        :raises ModelError: If one or more :class:`~gillespy2.core.events.Event` objects contain an invalid
            :attr:`Event.trigger` value. Specifically, it is asserted that :code:`Event.trigger` is not :code:`None` and
            the trigger contains the :code:`"expression"` property.

        :raises ParameterError: If some other error occurs while adding the :class:`~gillespy2.core.events.Event`
            object(s).
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
        Adds one or more :class:`~gillespy2.core.functiondefinition.FunctionDefinition` to this :class:`Model`
        instance.

        :param function_definitions: One or more function definition objects to be added to this model.

        :returns: One or more function definition objects that were successfully added to this :class:`Model`.
        :rtype: Union[FunctionDefinition, List[FunctionDefinition]]

        :raises ModelError: If :func:`Model.problem_with_name` invalidates the name of one or more
            :class:`~gillespy2.core.functiondefinition.FunctionDefinition` objects.

        :raises ParameterError:  If some other error occurs while adding the
            :class:`~gillespy2.core.functiondefinition.FunctionDefinition` object.
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
        Adds one or more :class:`~gillespy2.core.assignmentrule.AssignmentRule` objects to the :class:`Model`.

        :param assignment_rules: One or more :class:`~gillespy2.core.assignmentrule.AssignmentRule` objects. If multiple
            are provided then they must be within a list.

        :returns: One or more assignment rule objects that were successfully added to this :class:`Model`.
        :rtype: Union[AssignmentRule, List[AssignmentRule]]

        :raises ModelError: If the :code:`variable` property of a :class:`~gillespy2.core.assignmentrule.AssignmentRule` instance
            is already set within a :class:`~gillespy2.core.assignmentrule.AssignmentRule` or
            :class:`~gillespy2.core.raterule.RateRule` attached to this :class:`Model`.

        :raises ModelError: If one or more assignment rules in :code:`assignment_rules` has the same :code:`name` as a
            :class:`~gillespy2.core.assignmentrule.AssignmentRule` on this :class:`Model` instance.

        :raises ModelError: If one or more assignment rules in :code:`assignment_rules` has an invalid :attr:`AssignmentRule.formula`
            value. Specifically, if :code:`AssignmentRule.formula == ""`.

        :raises ModelError: If one or more assignment rules in :code:`assignment_rules` has a :attr:`AssignmentRule.variable` value of
            :code:`None`.

        :raises ModelError: If :func:`Model.problem_with_name` invalidates the name of one or more
            :class:`~gillespy2.core.assignmentrule.AssignmentRule` objects.

        :raises ModelError:
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

    def timespan(self, time_span: Union[TimeSpan, Iterable[int]]) -> None:
        """
        Sets the timespan of the simulation.

        .. note::
            StochKit **only** supports uniform timespans.

        :param time_span: The evently-spaced list of times at which :class:`~gillespy2.core.species.Species` populations
            are sampled during the simulation. It is recommended to use the :class:`~gillespy2.core timespan.TimeSpan`
            class instead of a list, :code:`numpy.ndarray`, or some other iterable collection.
        """
        if isinstance(time_span, TimeSpan) or type(time_span).__name__ == "TimeSpan":
            self.tspan = time_span
        else:
            self.tspan = TimeSpan(time_span)

    def get_event(self, ename: str) -> Event:
        """
        Returns an :class:`~gillespy2.core.events.Event` on this :class:`Model` by name.

        :param ename: The name of the event to return.

        :returns: An event with name :code:`ename`.
        :rtype: Event

        :raises KeyError: If a :class:`~gillespy2.core.events.Event` with name :code:`ename` does not exist on this
            :class:`Model`.
        """
        return self.listOfEvents[ename]

    def get_all_events(self) -> Dict[str, Event]:
        """
        Returns all :class:`~gillespy2.core.events.Event` objects attached to this :class:`Model` instance.

        :returns: A dictionary of events of the form :code:`{ name: Event }`.
        :rtype: Dict[str, Event]
        """
        return self.listOfEvents

    def delete_event(self, name: str) -> None:
        """
        Removes an :class:`~gillespy2.core.events.Event` object from this :class:`Model` by name.

        :param name: The name of the event object to remove.
        :raises KeyError: If an event with name :code:`name` does not exist.
        """
        self.listOfEvents.pop(name)
        if name in self._listOfEvents:
            self._listOfEvents.pop(name)

    def delete_all_events(self) -> None:
        """
        Removes all :class:`~gillespy2.core.events.Event` objects attached to this :class:`Model` instance.
        """
        self.listOfEvents.clear()
        self._listOfEvents.clear()

    def get_rate_rule(self, rname: str) -> RateRule:
        """
        Returns a :class:`~gillespy2.core.raterule.RateRule` on this :class:`Model` by name.

        :param p_name: The name of the rate rule to return.

        :returns: The rate rule object.
        :rtype: RateRule

        :raises KeyError: If a rate rule with name :class:`~gillespy2.core.raterule.RateRule` does not exist on this
            :class:`Model`.
        """
        return self.listOfRateRules[rname]

    def get_all_rate_rules(self) -> Dict[str, RateRule]:
        """
        Returns all :class:`~gillespy2.core.raterule.RateRule` objects attached to this :class:`Model`.

        :returns: A dictionary of event objects of the form :code:`{ name: RateRule }`.
        :rtype: Dict[str, RateRule]
        """
        return self.listOfRateRules

    def delete_rate_rule(self, name: str) -> None:
        """
        Removes a :class:`~gillespy2.core.raterule.RateRule` from this :class:`Model` by name.

        :param name: The name of the rate rule object to remove.
        :raises KeyError: If a rate rule with name :code:`name` does not exist.
        """
        self.listOfRateRules.pop(name)
        if name in self._listOfRateRules:
            self._listOfRateRules.pop(name)

    def delete_all_rate_rules(self) -> None:
        """
        Removes all :class:`~gillesp2.core.raterule.RateRule` objects attached to this :class:`Model` instance.
        """
        self.listOfRateRules.clear()
        self._listOfRateRules.clear()

    def get_assignment_rule(self, aname: str) -> AssignmentRule:
        """
        Returns an :class:`~gillespy2.core.assignmentrule.AssignmentRule` on this :class:`Model` by name.

        :param aname: The name of the assignment rule to return.

        :returns: An :class:`~gillespy2.core.assignmentrule.AssignmentRule` with name :code:`aname`.
        :rtype: AssignmentRule

        :raises KeyError: If an :class:`~gillespy2.core.assignmenrrule.AssignmentRule` with name :code:`aname` does not
            exist on this :class:`Model`.

        """
        return self.listOfAssignmentRules[aname]

    def get_all_assignment_rules(self) -> Dict[str, AssignmentRule]:
        """
        Returns all :class:`~gillespy2.core.assignmentrule.AssignmentRule` objects attached to this :class:`Model`
        instance.

        :returns A dictionary of assignment rules of the form :code:`{ name: AssignmentRule }`.
        :rtype: Dict[str, AssignmentRule]
        """
        return self.listOfAssignmentRules

    def delete_assignment_rule(self, name: str) -> None:
        """
        Removes an :class:`~gillespy2.core.assignmentrule.AssignmentRule` object from this :class:`Model` by name.

        :param name: The name of the assignment rule to remove.
        :raises KeyError: If an :class:`~gillespy2.core.assignmentrule.AssignmentRule` with name :code:`name` does not
            exist.
        """
        self.listOfAssignmentRules.pop(name)
        if name in self._listOfAssignmentRules:
            self._listOfAssignmentRules.pop(name)

    def delete_all_assignment_rules(self) -> None:
        """
        Removes all :class:`~gillespy2.core.assignmentrule.AssigmentRule` objects attached to this :class:`Model`.
        """
        self.listOfAssignmentRules.clear()
        self._listOfAssignmentRules.clear()

    def get_function_definition(self, fname: str) -> FunctionDefinition:
        """
        Returns a :class:`~gillespy2.functiondefinition.FunctionDefinition` on this :class:`Model` by name.

        :param fname: The name of the function definition to return.

        :returns: The function definition object.
        :rtype: FunctionDefinition

        :raises KeyError: If a :class:`~gillespy2.functiondefinition.FunctionDefinition` with name :code:`rname` does not exist on this :class:`Model`.
        """
        return self.listOfFunctionDefinitions[fname]

    def get_all_function_definitions(self) -> Dict[str, FunctionDefinition]:
        """
        Returns all :class:`~gillespy2.functiondefinition.FunctionDefinition` objects attached to this :class:`Model`.

        :returns: A dictionary of :class:`~gillespy2.functiondefinition.FunctionDefinition` objects of the form:
            :code:`{ name: FunctionDefinition }`.

        :rtype: Dict[str, FunctionDefinition]
        """
        return self.listOfFunctionDefinitions

    def delete_function_definition(self, name: str) -> None:
        """
        Removes a :class:`~gillespy2.functiondefinition.FunctionDefinition` object from this :class:`Model` by name.

        :param name: The name of the :class:`~gillespy2.functiondefinition.FunctionDefinition` to remove.
        :raises KeyError: If a :class:`~gillespy2.functiondefinition.FunctionDefinition` with name :code:`name` does
            not exist.
        """
        self.listOfFunctionDefinitions.pop(name)
        if name in self._listOfFunctionDefinitions:
            self._listOfFunctionDefinitions.pop(name)

    def delete_all_function_definitions(self) -> None:
        """
        Removes all :class:`~gillespy2.functiondefinition.FunctionDefinition` objects attached to this :class:`Model`
        instance.
        """
        self.listOfFunctionDefinitions.clear()
        self._listOfFunctionDefinitions.clear()

    def get_element(self, ename: str) -> Component:
        """
        Returns an element attached to this :class:`Model` instance by name.

        This function will return one of the following types:
            - :class:`gillespy2.Species`
            - :class:`gillespy2.Parameter`
            - :class:`gillespy2.Reaction`
            - :class:`gillespy2.Event`
            - :class:`gillespy2.RateRule`
            - :class:`gillespy2.AssignmentRule`
            - :class:`gillespy2.FunctionDefinition`
            - :class:`gillespy2.TimeSpan`

        :param ename: The name of the element to return.

        :returns: The element object.
        :rtype: One or more objects with type defined in the docstring above.

        :raises ModelError: If an element with name :code:`ename` does not exist on this :class:`Model`.
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
        Returns the best solver for the current :class:`Model` instance. Currently,
        :class:`~gillespy2.core.assignmentrule.AssignmentRule`, :class:`~gillespy2.core.raterule.RateRule`,
        :class:`~gillespy2.core.functiondefinition.FunctionDefinition`, :class:`~gillespy2.core.events.Event`, and
        :class:`~gillespy2.core.species.Species` with a dynamic or continuous population must use either
        :class:`~gillespy2.solvers.numpy.tau_hybrid_solver.TauHybridSolver` or
        :class:`~gillespy2.solvers.cpp.tau_hybrid_c_solver.TauHybridCSolver`.

        Optimized C++ solvers require :code:`gcc` and :code:`make` to be installed on the host system to run
        simulations. One can determine which dependencies are missing like so:

        .. code-block:: python

            from gillespy2.solvers.cpp.build.build_engine import BuildEngine

            print(BuildEngine.get_missing_dependencies())

        A model can still be simulated if dependencies are missing, however, with the :code:`numpy` suite of solvers.
        These are functionally identical to their C++ counterparts, but aren't as performant.

        The following cases can be used to determine which solver implementation will be chosen based on your Python
        environment.

        If :meth:`BuildEngine.get_missing_dependencies` returns a list of length :code:`0`:

            - :class:`~gillespy2.solvers.cpp.ode_c_solver.ODECSolver`
            - :class:`~gillespy2.solvers.cpp.ssa_c_solver.SSACSolver`
            - :class:`~gillespy2.solvers.cpp.tau_hybrid_c_solver.TauHybridCSolver`
            - :class:`~gillespy2.solvers.cpp.tau_leaping_c_solver.TauLeapingCSolver`

        If :meth:`BuildEngine.get_missing_dependencies` returns a list of length :code:`> 0`:

            - :class:`~gillespy2.solvers.numpy.CLE_solver.CLESolver`
            - :class:`~gillespy2.solvers.numpy.ode_solver.ODESolver`
            - :class:`~gillespy2.solvers.numpy.ssa_solver.NumPySSASolver`
            - :class:`~gillespy2.solvers.numpy.tau_hybrid_solver.TauHybridSolver`
            - :class:`~gillespy2.solvers.numpy.tau_leaping_solver.TauLeapingSolver`

        :returns: The chosen solver type.
        :rtype: A solver of inherited type :class:`~gillespy2.core.gillespySolver.GillesPySolver`.
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
        Returns the best solver implementation for the specified algorithm *type*.

        Specifically, :meth:`BuildEngine.get_missing_dependencies` is used to determine if the current Python
        environment (and by extension, the host machine) is missing C++ build dependencies. See
        :meth:`Model.get_best_solver` for a more thorough explanation on the usage of
        :class:`~gillespy2.solvers.cpp.build.build_engine.BuildEngine` and how to interpret its results.

        :param algorithm: The name of the algorithm to be used. Must be one of the following values:

            - :code:`"CLE"`
            - :code:`"ODE"`
            - :code:`"SSA"`
            - :code:`"Tau-Hybrid"`
            - :code:`"Tau-Leaping"`

        :returns: The solver type whose implementation is best suited for the current Python environment.
        :rtype: Refer to :meth:`Model.get_best_solver` for a thorough overview of what may be returned.

        :raises ModelError: If neither :code:`numpy` or C++ build dependencies are installed.
        :raises ModelError: If :code:`algorithm` is of some other value.
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
        Returns a list of solver-specific features that are present within this :class:`Model`. This is used to
        determine if the current :class:`Model` is compatible with a solver.

        :returns: A set containing the classes of every solver-specific feature present within this model.
        :rtype: AbstractSet[Union[Event, RateRule, AssignmentRule, FunctionDefinition]]
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
        Creates an StochKit XML document from an exisiting :class:`Model` object.
        This method assumes that all the parameters in the model are already
        resolved to scalar floats (see Model.resolveParamters).

        Note, this method is intended to be used internally by the models
        'serialization' function, which performs additional operations and
        tests on the model prior to writing out the XML file.

        You should NOT do:

        .. code-block:: python

            document = StochMLDocument.fromModel(model)
            print(document.toString())

        You SHOULD do:

        .. code-block:: python

            print(model.serialize())

        :param model: The :class:`Model` instance to convert to a StochML document.

        :returns: The newly created StochML document.
        :rtype: StochMLDocument

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
        """
        Reads, parses, and returns a StochML document from a file on disk.

        :param filepath: The path of the StochML document on disk.

        :returns: The parsed StochML document.
        :rtype: StochMLDocument
        """
        tree = eTree.parse(filepath)
        root = tree.getroot()
        md = cls()
        md.document = root
        return md

    @classmethod
    def from_string(cls, string: str) -> "StochMLDocument":
        """
        Returns a StochML document parsed from the specified string's contents.

        :param string: The string to be parsed.

        :returns: The parsed StochMLDocument.
        :rtype: StochMLDocument
        """
        root = eTree.fromString(string)

        md = cls()
        md.document = root
        return md

    def to_model(self, name: str) -> Model:
        """
        Returns the gillespy2 :class:`Model` representation of the current :class:`StochMLDocument` instance.

        :param name: The name of the newly created :class:`Model` object. If unset, the current :class:`StochMLDocument`
            instance will be parsed for a name instead.

        :returns: The newly created :class:`Model` object.
        :rtype: Model
        """

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
        """
        Returns the string representation of the current :class:`StochMLDocument` instance.

        :returns: The string representation of the current StochMLDocument.
        :rtype: str
        """
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
