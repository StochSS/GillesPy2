from gillespy2.core.reaction import *
from gillespy2.core.raterule import RateRule
from gillespy2.core.parameter import Parameter
from gillespy2.core.species import Species
from gillespy2.core.reaction import Reaction
import numpy as np
from gillespy2.core.results import Trajectory,Results
from collections import OrderedDict
from gillespy2.core.gillespyError import *

try:
    import lxml.etree as eTree

    no_pretty_print = False

except:
    import xml.etree.ElementTree as eTree
    import xml.dom.minidom
    import re
    no_pretty_print = True


def import_SBML(filename, name=None, gillespy_model=None):
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
    """

    try:
        from gillespy2.sbml.SBMLimport import convert
    except ImportError:
        raise ImportError('SBML conversion not imported successfully')

    return convert(filename, model_name=name, gillespy_model=gillespy_model)


class Model(SortableObject):
    # reserved names for model species/parameter names, volume, and operators.
    reserved_names = ['vol']
    special_characters = ['[', ']', '+', '-', '*', '/', '.', '^']

    """
    Representation of a well mixed biochemical model. Contains reactions,
    parameters, species.

    :param name: The name of the model, or an annotation describing it.
    :type name: str
    
    :param population: The type of model being described. A discrete stochastic model is a
    population model (True), a deterministic model is a concentration model
    (False). Automatic conversion from population to concentration models
    may be used, by setting the volume parameter.
    Attributes
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
                raise Warning(
                    "Concentration models account for volume implicitly, explicit volume definition is not required. "
                    "Note: concentration models may only be simulated deterministically.")

        self.volume = volume

        # Dict that holds flattended parameters and species for
        # evaluation of expressions in the scope of the model.
        self.namespace = OrderedDict([])

        if tspan is None:
            self.timespan(np.linspace(0, 20, 401))
        else:
            self.timespan(tspan)

    def __str__(self):
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
        return print_string

    def serialize(self):
        """ Serializes the Model object to valid StochML. """
        self.resolve_parameters()
        doc = StochMLDocument().from_model(self)
        return doc.to_string()

    def update_namespace(self):
        """ Create a dict with flattened parameter and species objects. """
        self.namespace = OrderedDict([])
        for param in self.listOfParameters:
            self.namespace[param] = self.listOfParameters[param].value

    def sanitized_species_names(self):
        """
        Generate a dictionary mapping user chosen species names to simplified formats which will be used
        later on by GillesPySolvers evaluating reaction propensity functions.
        :return: the dictionary mapping user species names to their internal GillesPy notation.
        """
        species_name_mapping = OrderedDict([])
        for i, name in enumerate(sorted(self.listOfSpecies.keys())):
            species_name_mapping[name] = 'S[{}]'.format(i)
        return species_name_mapping

    def problem_with_name(self, name):
        if name in Model.reserved_names:
            return ModelError(
                'Name "{}" is unavailable. It is reserved for internal GillesPy use. Reserved Names: ({}).'.format(name,
                                                                                                                   Model.reserved_names))
        if name in self.listOfSpecies:
            return ModelError('Name "{}" is unavailable. A species with that name exists.'.format(name))
        if name in self.listOfParameters:
            return ModelError('Name "{}" is unavailable. A parameter with that name exists.'.format(name))
        if name.isdigit():
            return ModelError('Name "{}" is unavailable. Names must not be numeric strings.'.format(name))
        for special_character in Model.special_characters:
            if special_character in name:
                return ModelError(
                    'Name "{}" is unavailable. Names must not contain special characters: {}.'.format(name,
                                                                                                      Model.special_characters))

    def get_species(self, s_name):
        """
        Returns a species object by name.

        :param s_name: Name of the species object to be returned:
        :type s_name: str
        """
        return self.listOfSpecies[s_name]

    def get_all_species(self):
        """
        Returns a dict of all species in the model, of the form:
        {name : species object}
        """
        return self.listOfSpecies

    def add_species(self, obj):
        """
        Adds a species, or list of species to the model.

        :param obj: The species or list of species to be added to the model object
        :type obj: Species, or list of species
        """

        if isinstance(obj, list):
            for S in sorted(obj):
                self.add_species(S)
        else:
            try:
                problem = self.problem_with_name(obj.name)
                if problem is not None:
                    raise problem
                self.listOfSpecies[obj.name] = obj
                self._listOfSpecies[obj.name] = 'S{}'.format(len(self._listOfSpecies))
            except Exception as e:
                raise ParameterError("Error using {} as a Species. Reason given: {}".format(obj, e))
        return obj

    def delete_species(self, obj):
        """
        Removes a species object by name.

        :param obj: Name of the species object to be removed
        :type obj: str
        """
        self.listOfSpecies.pop(obj)
        self._listOfSpecies.pop(obj)

    def delete_all_species(self):
        """
        Removes all species from the model object.
        """
        self.listOfSpecies.clear()
        self._listOfSpecies.clear()

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

    def sanitized_parameter_names(self):
        """
        Generate a dictionary mapping user chosen parameter names to simplified formats which will be used
        later on by GillesPySolvers evaluating reaction propensity functions.
        :return: the dictionary mapping user parameter names to their internal GillesPy notation.
        """
        parameter_name_mapping = OrderedDict()
        parameter_name_mapping['vol'] = 'V'
        for i, name in enumerate(sorted(self.listOfParameters.keys())):
            if name not in parameter_name_mapping:
                parameter_name_mapping[name] = 'P{}'.format(i)
        return parameter_name_mapping

    def get_parameter(self, p_name):
        """
        Returns a parameter object by name.

        :param p_name: Name of the parameter object to be returned
        :type p_name: str
        """
        try:
            return self.listOfParameters[p_name]
        except:
            raise ModelError("No parameter named " + p_name)

    def get_all_parameters(self):
        """
        Returns a dict of all parameters in the model, of the form:
        {name : parameter object}
        """
        return self.listOfParameters

    def add_parameter(self, params):
        """
        Adds a parameter, or list of parameters to the model.

        :param params:  The parameter or list of parameters to be added to the model object.
        :type params: Parameter, or list of parameters
        """
        if isinstance(params, list):
            for p in sorted(params):
                self.add_parameter(p)
        else:
            try:
                problem = self.problem_with_name(params.name)
                if problem is not None:
                    raise problem
                self.listOfParameters[params.name] = params
                self._listOfParameters[params.name] = 'P{}'.format(len(self._listOfParameters))
            except Exception as e:
                raise ParameterError("Error using {} as a Parameter. Reason given: {}".format(params, e))
        return params

    def delete_parameter(self, obj):
        """
        Removes a parameter object by name.

        :param obj: Name of the parameter object to be removed
        :type obj: str
        """
        self.listOfParameters.pop(obj)
        self._listOfParameters.pop(obj)

    def set_parameter(self, p_name, expression):
        """
        Set the value of an existing parameter "pname" to "expression".

        :param p_name: Name of the parameter whose value will be set.
        :type p_name: str

        :param expression: String that may be executed in C, describing the value of the
        parameter. May reference other parameters by name. (e.g. "k1*4")
        :type expression: str
        """

        p = self.listOfParameters[p_name]
        p.expression = expression
        p.evaluate()

    def resolve_parameters(self):
        """ Internal function:
        attempt to resolve all parameter expressions to scalar floats.
        This methods must be called before exporting the model.
        """
        self.update_namespace()
        for param in self.listOfParameters:
            try:
                self.listOfParameters[param].evaluate(self.namespace)
            except:
                raise ParameterError("Could not resolve Parameter expression {} to a scalar value.".format(param))

    def delete_all_parameters(self):
        """ Deletes all parameters from model. """
        self.listOfParameters.clear()
        self._listOfParameters.clear()

    def validate_reactants_and_products(self, reactions):
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

    def add_reaction(self, reactions):
        """
        Adds a reaction, or list of reactions to the model.

        :param reactions: The reaction or list of reaction objects to be added to the model
        object.
        :type reactions: Reaction, or list of Reactions
        """

        # TODO, make sure that you cannot overwrite an existing reaction
        if isinstance(reactions, list):
            for r in sorted(reactions):
                self.add_reaction(r)
        else:
            try:
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

    def add_rate_rule(self, rate_rules):
        """
        Adds a rate rule, or list of rate rules to the model.

        :param rate_rules: The rate rule or list of rate rule objects to be added to the model object.
        :type rate_rules: RateRule, or list of RateRules
        """
        if isinstance(rate_rules, list):
            for rr in sorted(rate_rules):
                self.add_rate_rule(rr)
        else:
            try:
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

                self.listOfRateRules[rate_rules.name] = rate_rules
                sanitized_rate_rule = RateRule(name='RR{}'.format(len(self._listOfRateRules)))
                sanitized_rate_rule.formula = rate_rules.sanitized_formula(self._listOfSpecies,
                                                                           self._listOfParameters)
                self._listOfRateRules[rate_rules.name] = sanitized_rate_rule
            except Exception as e:
                raise ParameterError("Error using {} as a Rate Rule. Reason given: {}".format(rate_rules, e))
        return rate_rules

    def add_event(self, event):
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

    def add_function_definition(self, function_definitions):
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
                self.listOfFunctionDefinitions[function_definitions.name] = function_definitions
            except Exception as e:
                raise ParameterError(
                    "Error using {} as a Function Definition. Reason given: ".format(function_definitions, e))

    def add_assignment_rule(self, assignment_rules):
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
                raise ParameterError("Error using {} as a Assignment Rule. Reason given: ".format(assignment_rules, e))

    def timespan(self, time_span):
        """
        Set the time span of simulation. StochKit does not support non-uniform
        timespans.

        :param time_span: Evenly-spaced list of times at which to sample the species
        populations during the simulation.
        :type time_span: numpy ndarray
        """

        items = np.diff(time_span)
        items = np.array([round(item, 10) for item in items])
        isuniform = (len(set(items)) == 1)

        if isuniform:
            self.tspan = time_span
        else:
            raise InvalidModelError("StochKit only supports uniform timespans")

    def get_reaction(self, rname):
        """
        :param rname: name of reaction to return
        :return: Reaction object
        """
        return self.listOfReactions[rname]

    def get_all_reactions(self):
        """
        :return: dict of all Reaction objects
        """
        return self.listOfReactions

    def delete_reaction(self, obj):
        """
        :param obj: Name of Reaction to be removed
        """
        self.listOfReactions.pop(obj)
        self._listOfReactions.pop(obj)

    def delete_all_reactions(self):
        """
        Clears all reactions in model
        """
        self.listOfReactions.clear()
        self._listOfReactions.clear()

    def get_event(self, ename):
        """
        :param ename: Name of Event to get
        :return: Event object
        """
        return self.listOfEvents[ename]

    def get_all_events(self):
        """
        :return: dict of all Event objects
        """
        return self.listOfEvents

    def delete_event(self, ename):
        """
        Removes specified Event from model
        :param ename: Name of Event to be removed
        """
        self.listOfEvents.pop(ename)
        self._listOfEvents.pop(ename)

    def delete_all_events(self):
        """
        Clears models events
        """
        self.listOfEvents.clear()
        self._listOfEvents.clear()

    def get_rate_rule(self, rname):
        """
        :param rname: Name of Rate Rule to get
        :return: RateRule object
        """
        return self.listOfRateRules[rname]

    def get_all_rate_rules(self):
        """
        :return: dict of all Rate Rule objects
        """
        return self.listOfRateRules

    def delete_rate_rule(self, rname):
        """
        Removes specified Rate Rule from model
        :param rname: Name of Rate Rule to be removed
        """
        self.listOfRateRules.pop(rname)
        self._listOfRateRules.pop(rname)

    def delete_all_rate_rules(self):
        """
        Clears all of models Rate Rules
        """
        self.listOfRateRules.clear()
        self._listOfRateRules.clear()

    def get_assignment_rule(self, aname):
        """
        :param aname: Name of Assignment Rule to get
        :return: Assignment Rule object
        """
        return self.listOfAssignmentRules[aname]

    def get_all_assignment_rules(self):
        """
        :return: dict of models Assignment Rules
        """
        return self.listOfAssignmentRules

    def delete_assignment_rule(self, aname):
        """
        Removes an assignment rule from a model
        :param aname: Name of AssignmentRule object to be removed from model
        """
        self.listOfAssignmentRules.pop(aname)
        self._listOfAssignmentRules.pop(aname)

    def delete_all_assignment_rules(self):
        """
        Clears all assignment rules from model
        """
        self.listOfAssignmentRules.clear()
        self._listOfAssignmentRules.clear()

    def get_function_definition(self, fname):
        """
        :param fname: name of Function to get
        :return: FunctionDefinition object
        """
        return self.listOfFunctionDefinitions[fname]

    def get_all_function_definitions(self):
        """
        :return: Dict of models function definitions
        """
        return self.listOfFunctionDefinitions

    def delete_function_definition(self, fname):
        """
        Removes specified Function Definition from model
        :param fname: Name of Function Definition to be removed
        """
        self.listOfFunctionDefinitions.pop(fname)
        self._listOfFunctionDefinitions.pop(fname)

    def delete_all_function_definitions(self):
        """
        Clears all Function Definitions from a model
        """
        self.listOfFunctionDefinitions.clear()
        self._listOfFunctionDefinitions.clear()

    def get_element(self, ename):
        """
        get element specified by name
        :param ename: name of element to search for
        :return:value of element, or 'element not found'
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
        return 'Element not found!'


    def get_best_solver(self, precompile=True):
        """
        Finds best solver for the users simulation. Currently, AssignmentRules, RateRules, FunctionDefinitions,
        Events, and Species with a dynamic, or continuous population must use the TauHybridSolver.
        :param precompile: If True, and the model contains no AssignmentRules, RateRules, FunctionDefinitions, Events,
        or Species with a dynamic or continuous population, the get_best_solver will choose the VariableSSACSolver, else
        it will choose SSACSolver
        :type precompile: bool
        :return: gillespy2.gillespySolver
        """
        from gillespy2.solvers.numpy import can_use_numpy
        hybrid_check = False
        if len(self.get_all_assignment_rules()) or len(self.get_all_rate_rules())  \
                or len(self.get_all_function_definitions()) or len(self.get_all_events()):
            hybrid_check = True

        if len(self.get_all_species()) and hybrid_check == False:
            for i in self.get_all_species():
                tempMode = self.get_species(i).mode
                if tempMode == 'dynamic' or tempMode == 'continuous':
                    hybrid_check = True
                    break
        if can_use_numpy and hybrid_check:
            from gillespy2 import TauHybridSolver
            return TauHybridSolver

        elif not can_use_numpy and hybrid_check:
            raise ModelError('TauHybridSolver is the only solver currently that supports '
                             'AssignmentRules, RateRules, FunctionDefinitions, or Events. '
                             'Please install Numpy.')
        from gillespy2.solvers.utilities.cpp_support_test import cpp_support

        if cpp_support is False and can_use_numpy and not hybrid_check:
            from gillespy2 import NumPySSASolver
            return NumPySSASolver

        else:
            if precompile:
                from gillespy2 import VariableSSACSolver
                return VariableSSACSolver
            else:
                from gillespy2 import SSACSolver
                return SSACSolver

    def run(self, solver=None, timeout=0, t=None, show_labels=True, cpp_support=False, **solver_args):
        """
        Function calling simulation of the model. There are a number of
        parameters to be set here.

        :param solver: The solver by which to simulate the model. This solver object may
        be initialized separately to specify an algorithm. Optional,
        defaults to ssa solver.
        :type solver: gillespy.GillesPySolver

        :param timeout: Allows a time_out value in seconds to be sent to a signal handler, restricting simulation run-time
        :type timeout: int

        :param t: End time of simulation
        :type t: int
        :param solver_args: Solver-specific arguments to be passed to solver.run()

        :param cpp_support: INTERNAL USE ONLY, flag for whether or not a computer has the capability to compile a
        C++ program.
        :type cpp_support: bool

        :return  If show_labels is False, returns a numpy array of arrays of species population data. If show_labels is
        True,returns a Results object that inherits UserList and contains one or more Trajectory objects that
        inherit UserDict. Results object supports graphing and csv export.

        To pause a simulation and retrieve data before the simulation, keyboard interrupt the simulation by pressing
        control+c or pressing stop on a jupyter notebook. To resume a simulation, pass your previously ran results
        into the run method, and set t = to the time you wish the resuming simulation to end (run(resume=results, t=x)).
        Pause/Resume is only supported for SINGLE TRAJECTORY simulations. T MUST BE SET OR UNEXPECTED BEHAVIOR MAY OCCUR.
        """

        if not show_labels:
            from gillespy2.core import log
            log.warning('show_labels = False is deprecated. Future releases '
                        'of GillesPy2 may not support this feature.')
        if t is None:
            t = self.tspan[-1]

        if solver is None:
            solver = self.get_best_solver()

        try:
            solver_results, rc = solver.run(model=self, t=t, increment=self.tspan[-1] - self.tspan[-2],
                                            timeout=timeout, **solver_args)
        except Exception as e:
            # If user has specified the SSACSolver, but they don't actually have a g++ compiler,
            # This will throw an error and throw log. IF a user specifies cpp_support == True and don't have a compiler
            # They would bypass this log.warning and just recieve an error
            if cpp_support is False and not isinstance(solver, str):
                if solver.name == 'SSACSolver' or solver.name == 'VariableSSACSolver':
                    from gillespy2.core import log
                    log.warning("Please install/configure 'g++' and 'make' on your"
                                " system, to ensure that GillesPy2 C solvers will"
                                " run properly.")
            raise SimulationError(
                "argument 'solver={}' to run() failed.  Reason Given: {}".format(solver, e))

        if rc == 33:
            from gillespy2.core import log
            log.warning('GillesPy2 simulation exceeded timeout.')

        if hasattr(solver_results[0], 'shape'):
            return solver_results

        if len(solver_results) > 0:
            results_list = []
            for i in range(0, len(solver_results)):
                temp = Trajectory(data=solver_results[i], model=self, solver_name=solver.name, rc=rc)
                results_list.append(temp)

            results = Results(results_list)
            if show_labels == False:
                results = results.to_array()
            return results

        if len(solver_results) > 0:
            results_list = []
            for i in range(0, len(solver_results)):
                temp = Trajectory(data=solver_results[i], model=self, solver_name=solver.name, rc=rc)
                results_list.append(temp)

            results = Results(results_list)
            if show_labels == False:
                results = results.to_array()
            return results

        else:
            raise ValueError("number_of_trajectories must be non-negative and non-zero")


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
        Creates an StochKit XML document from an exisiting Mdoel object.
        This method assumes that all the parameters in the model are already
        resolved to scalar floats (see Model.resolveParamters).

        Note, this method is intended to be used interanally by the models
        'serialization' function, which performs additional operations and
        tests on the model prior to writing out the XML file. You should NOT \
        do:

        document = StochMLDocument.fromModel(model)
        print document.toString()

        You SHOULD do

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
    def from_file(cls, filepath):
        """ Intializes the document from an exisiting native StochKit XML
        file read from disk. """
        tree = eTree.parse(filepath)
        root = tree.getroot()
        md = cls()
        md.document = root
        return md

    @classmethod
    def from_string(cls, string):
        """ Intializes the document from an exisiting native StochKit XML
        file read from disk. """
        root = eTree.fromString(string)

        md = cls()
        md.document = root
        return md

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
                p.evaluate()
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
                        p.evaluate()
                        model.add_parameter(p)
                        reaction.marate = model.listOfParameters[
                            generated_rate_name]

                    reaction.__create_mass_action()
                except Exception as e:
                    raise
            elif type == 'customized':
                try:
                    propfunc = reac.find('PropensityFunction').text
                except Exception as e:
                    raise InvalidStochMLError(
                        "Found a customized propensity function, but no expression was given. {}".format(e))
                reaction.propensity_function = propfunc
            else:
                raise InvalidStochMLError(
                    "Unsupported or no reaction type given for reaction" + name)

            model.add_reaction(reaction)

        return model

    def to_string(self):
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

    def __species_to_element(self, S):
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

    def __parameter_to_element(self, P):
        e = eTree.Element('Parameter')
        idElement = eTree.Element('Id')
        idElement.text = P.name
        e.append(idElement)
        expressionElement = eTree.Element('Expression')
        expressionElement.text = str(P.value)
        e.append(expressionElement)
        return e

    def __reaction_to_element(self, R, model_volume):
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
            srElement.set('id', str(reactant))
            srElement.set('stoichiometry', str(stoichiometry))
            reactants.append(srElement)

        e.append(reactants)

        products = eTree.Element('Products')
        for product, stoichiometry in R.products.items():
            srElement = eTree.Element('SpeciesReference')
            srElement.set('id', str(product))
            srElement.set('stoichiometry', str(stoichiometry))
            products.append(srElement)
        e.append(products)

        return e
