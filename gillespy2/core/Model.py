from gillespy2.core.SortableObject import SortableObject
from gillespy2.core.gillespyError import *
from gillespy2.core.Reaction import *
from gillespy2.core.RateRule import RateRule
from gillespy2.core.results import Trajectory,Results
from gillespy2.core.Utilities import StochMLDocument
from collections import OrderedDict
import numpy as np


class Model(SortableObject):
    # reserved names for model species/parameter names, volume, and operators.
    reserved_names = ['vol']
    special_characters = ['[', ']', '+', '-', '*', '/', '.', '^']

    """
    Representation of a well mixed biochemical model. Contains reactions,
    parameters, species.

    Attributes
    ----------
    name : str
        The name of the model, or an annotation describing it.
    population : bool
        The type of model being described. A discrete stochastic model is a
        population model (True), a deterministic model is a concentration model
        (False). Automatic conversion from population to concentration models
        may be used, by setting the volume parameter.
    volume : float
        The volume of the system matters when converting to from population to
        concentration form. This will also set a parameter "vol" for use in
        custom (i.e. non-mass-action) propensity functions.
    tspan : numpy ndarray
        The timepoints at which the model should be simulated. If None, a
        default timespan is added. May be set later, see Model.timespan
    annotation : str (optional)
        Optional further description of model
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

        Attributes
        ----------
        s_name : str
            Name of the species object to be returned.
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

        Attributes
        ----------
        obj : Species, or list of Species
            The species or list of species to be added to the model object.
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

        Attributes
        ----------
        obj : str
            Name of the species object to be removed.
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

        Attributes
        ----------
        units : str
            Either "population" or "concentration"
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

        Attributes
        ----------
        p_name : str
            Name of the parameter object to be returned.
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

        Attributes
        ----------
        obj : Parameter, or list of Parameters
            The parameter or list of parameters to be added to the model object.
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

        Attributes
        ----------
        obj : str
            Name of the parameter object to be removed.
        """
        self.listOfParameters.pop(obj)
        self._listOfParameters.pop(obj)

    def set_parameter(self, p_name, expression):
        """
        Set the value of an existing paramter "pname" to "expression".

        Attributes
        ----------
        p_name : str
            Name of the parameter whose value will be set.
        expression : str
            *String* that may be executed in C, describing the value of the
            parameter. May reference other parameters by name. (e.g. "k1*4")
        """

        p = self.listOfParameters[p_name]
        p.expression = expression
        p.evaluate()

    def resolve_parameters(self):
        """ Internal function:
        attempt to resolve all parameter expressions to scalar floats.
        This methods must be called before exporting the model. """
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
        for reactant in reactions.reactants.keys():
            if isinstance(reactant, str):
                if reactant not in self.listOfSpecies.keys():
                    raise ModelError(
                        'reactant: {0} for reaction {1} -- not found in model.listOfSpecies'.format(reactant,
                                                                                                    reactions.name))
                reactions.reactants[self.listOfSpecies[reactant]] = reactions.reactants[reactant]
                del reactions.reactants[reactant]
        for product in reactions.products.keys():
            if isinstance(product, str):
                if product not in self.listOfSpecies.keys():
                    raise ModelError('product: {0} for reaction {1} -- not found in model.listOfSpecies'.format(product,
                                                                                                                reactions.name))
                reactions.products[self.listOfSpecies[product]] = reactions.products[product]
                del reactions.products[product]

    def add_reaction(self, reactions):
        """
        Adds a reaction, or list of reactions to the model.

        Attributes
        ----------
        obj : Reaction, or list of Reactions
            The reaction or list of reaction objects to be added to the model
            object.
        """

        # TODO, make sure that you cannot overwrite an existing reaction
        if isinstance(reactions, list):
            for r in sorted(reactions):
                self.add_reaction(r)
        else:
            try:
                reactions.verify()
                self.validate_reactants_and_products(reactions)
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

                Attributes
                ----------
                obj : RateRule, or list of RateRules
                    The rate rule or list of rate rule objects to be added to the model
                    object.
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

                Attributes
                ----------
                event : Event, or list of Events
                    The event or list of event objects to be added to the model
                    object.
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
                        if a.variable in self.listOfSpecies:
                            a.variable = self.listOfSpecies[a.variable]
                        else:
                            raise ModelError('{0} not a valid Species'.format(a.variable))
                self.listOfEvents[event.name] = event
            except Exception as e:
                raise ParameterError("Error using {} as Event. Reason given: {}".format(event, e))
        return event

    def add_function_definition(self, function_definitions):
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

        tspan : numpy ndarray
            Evenly-spaced list of times at which to sample the species
            populations during the simulation.
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

    def run(self, solver=None, timeout=0, **solver_args):
        """
        Function calling simulation of the model. There are a number of
        parameters to be set here.

        Return
        ----------

        If show_labels is False, returns a numpy array of arrays of species population data. If show_labels is
        True,returns a Results object that inherits UserList and contains one or more Trajectory objects that
        inherit UserDict. Results object supports graphing and csv export.

        Attributes
        ----------
        solver : gillespy.GillesPySolver
            The solver by which to simulate the model. This solver object may
            be initialized separately to specify an algorithm. Optional,
            defaults to ssa solver.
        timeout : int
            Allows a time_out value in seconds to be sent to a signal handler, restricting simulation run-time
        solver_args :
            solver-specific arguments to be passed to solver.run()
        """

        if solver is not None:
            try:
                solver_results, rc = solver.run(model=self, t=self.tspan[-1],
                                                increment=self.tspan[-1] - self.tspan[-2], timeout=timeout,
                                                **solver_args)
            except Exception as e:
                raise SimulationError(
                    "argument 'solver={}' to run() failed.  Reason Given: {}".format(solver, e))
        else:
            from gillespy2.solvers.auto import SSASolver
            solver = SSASolver
            solver_results, rc = SSASolver.run(model=self, t=self.tspan[-1],
                                               increment=self.tspan[-1] -
                                                         self.tspan[-2], timeout=timeout, **solver_args)

        if rc == 33:
            from gillespy2.core import log
            log.warning('GillesPy2 simulation exceeded timeout.')

        if hasattr(solver_results[0], 'shape'):
            return solver_results
        if len(solver_results) is 1:
            results_list = []
            results_list.append(Trajectory(data=solver_results[0], model=self,
                                           solver_name=solver.name, rc=rc))
            return Results(results_list)

        if len(solver_results) > 1:
            results_list = []
            for i in range(0, solver_args.get('number_of_trajectories')):
                results_list.append(Trajectory(data=solver_results[i], model=self, solver_name=solver.name,
                                               rc=rc))
            return Results(results_list)


        else:
            raise ValueError("number_of_trajectories must be non-negative and non-zero")