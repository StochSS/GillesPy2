"""
A simple toolkit for creating and simulating discrete stochastic models in
python.

"""
from __future__ import division
from collections import OrderedDict
from gillespy2.core.results import Results,EnsembleResults
from gillespy2.core.gillespySolver import GillesPySolver
from gillespy2.core.gillespyError import *
import numpy as np

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

    Attributes
    ----------
    filename : str
        Path to the SBML file for conversion.
    name : str
        Name of the resulting model.
    gillespy_model : gillespy.Model
        If desired, the SBML model may be added to an existing GillesPy model.
    """

    try:
        from gillespy2.sbml.SBMLimport import convert
    except ImportError:
        raise ImportError('SBML conversion not imported successfully')

    return convert(filename, model_name=name, gillespy_model=gillespy_model)


class SortableObject(object):
    """Base class for GillesPy2 objects that are sortable."""

    def __eq__(self, other):
        return (isinstance(other, self.__class__)
                and ordered(self) == ordered(other))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        return not __le__(self, other)

    def __ge__(self, other):
        return not __lt__(self, other)

    def __lt__(self, other):
        if hasattr(self, 'id') and hasattr(other, 'id'):
            return self.id.lower() < other.id.lower()
        elif hasattr(self, 'name') and hasattr(other, 'name'):
            return self.name.lower() < other.name.lower()
        else:
            return repr(self) < repr(other)

    def __le__(self, other):
        if hasattr(self, 'id') and hasattr(other, 'id'):
            return self.id.lower() <= other.id.lower()
        elif hasattr(self, 'name') and hasattr(other, 'name'):
            return self.name.lower() <= other.name.lower()
        else:
            return repr(self) <= repr(other)

    def __cmp__(self, other):
        if hasattr(self, 'id') and hasattr(other, 'id'):
            return cmp(self.id.lower(), other.id.lower())
        elif hasattr(self, 'name') and hasattr(other, 'name'):
            return cmp(self.name.lower(), other.name.lower())
        else:
            return cmp(repr(self), repr(other))

    def __hash__(self):
        if hasattr(self, '_hash'):
            return self._hash
        if hasattr(self, 'id'):
            self._hash = hash(self.id)
        elif hasattr(self, 'name'):
            self._hash = hash(self.name)
        else:
            self._hash = hash(self)
        return self._hash


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

        # Dictionaries with Species, Reactions and Parameter objects.
        # Species, Reaction and Parameter names are used as keys.
        self.listOfParameters = OrderedDict()
        self.listOfSpecies = OrderedDict()
        self.listOfReactions = OrderedDict()
        self.listOfRateRules = OrderedDict()

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
            return ModelError('Name "{}" is unavailable. It is reserved for internal GillesPy use. Reserved Names: ({}).'.format(name, Model.reserved_names))
        if name in self.listOfSpecies:
            return ModelError('Name "{}" is unavailable. A species with that name exists.'.format(name))
        if name in self.listOfParameters:
            return ModelError('Name "{}" is unavailable. A parameter with that name exists.'.format(name))
        if name.isdigit():
            return ModelError('Name "{}" is unavailable. Names must not be numeric strings.'.format(name))
        for special_character in Model.special_characters:
            if special_character in name:
                return ModelError('Name "{}" is unavailable. Names must not contain special characters: {}.'.format(name, Model.special_characters))

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

        if isinstance(obj, Species):
            problem = self.problem_with_name(obj.name)
            if problem is not None:
                raise problem
            self.listOfSpecies[obj.name] = obj
        elif isinstance(obj, list):
            for S in sorted(obj):
                self.add_species(S)
        else:
            raise ModelError("Unexpected parameter for add_species. Parameter must be Species or list of Species.")
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

    def delete_all_species(self):
        """
        Removes all species from the model object.
        """
        self.listOfSpecies.clear()

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
        if isinstance(params,list):
            for p in sorted(params):
                self.add_parameter(p)
        else:
            if isinstance(params, Parameter):
                problem = self.problem_with_name(params.name)
                if problem is not None:
                    raise problem
                self.listOfParameters[params.name] = params
            else:
                raise ParameterError("Could not resolve Parameter expression {} to a scalar value.".format(params))
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

    def validate_reactants_and_products(self, reactions):
            for reactant in reactions.reactants.keys():
                if isinstance(reactant, str):
                    if reactant not in self.listOfSpecies.keys():
                        raise ModelError('reactant: {0} for reaction {1} -- not found in model.listOfSpecies'.format(reactant, reactions.name))
                    reactions.reactants[self.listOfSpecies[reactant]] = reactions.reactants[reactant]
                    del reactions.reactants[reactant]
            for product in reactions.products.keys():
                if isinstance(product, str):
                    if product not in self.listOfSpecies.keys():
                        raise ModelError('product: {0} for reaction {1} -- not found in model.listOfSpecies'.format(product, reactions.name))
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
        if isinstance(reactions,list):
            for r in sorted(reactions):
                self.add_reaction(r)
        elif isinstance(reactions,Reaction):
            reactions.verify()
            self.validate_reactants_and_products(reactions)
            if reactions.name in self.listOfReactions:
                raise ModelError("Duplicate name of reaction: {0}".format(reactions.name))
            self.listOfReactions[reactions.name] = reactions
        else:
            raise ModelError("Unexpected parameter for add_reaction. Parameter must be Reaction or list of Reactions.")
        return reactions

    def add_rate_rule(self, rate_rules):
        """
                Adds a rate rule, or list of rate rules to the model.

                Attributes
                ----------
                obj : RateRule, or list of RateRules
                    The reaction or list of raterule objects to be added to the model
                    object.
                """

        # TODO, make sure that you cannot overwrite an existing reaction
        # param_type = type(reactions).__name__
        if isinstance(rate_rules, list):
            for rr in sorted(rate_rules):
                self.add_rate_rule(rr)
        elif isinstance(rate_rules, RateRule):
            if rate_rules.species is None or not isinstance(rate_rules.species, Species): raise ModelError(
                'A Rate Rule must be associated with a valid species.')
            if rate_rules.expression == '': raise ModelError('Invalid Rate Rule. Expression must be a non-empty string value')
            self.listOfRateRules[rate_rules.species.name] = rate_rules
        else:
            raise ParameterError("Add_rate_rule accepts a RateRule object or a List of RateRule Objects")
        return rate_rules

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
        return self.listOfReactions[rname]

    def get_all_reactions(self):
        return self.listOfReactions

    def delete_reaction(self, obj):
        self.listOfReactions.pop(obj)

    def delete_all_reactions(self):
        self.listOfReactions.clear()

    def run(self, solver=None, **solver_args):
        """
        Function calling simulation of the model. There are a number of
        parameters to be set here.

        Return
        ----------

        If show_labels is False, returns a numpy array of arrays of species population data. If show_labels is True and
        number_of_trajectories is 1, returns a results object that inherits UserDict and supports plotting functions.
        If show_labels is False and number_of_trajectories is greater than 1, returns an ensemble_results object that
        inherits UserList and contains results objects and supports ensemble graphing.

        Attributes
        ----------
        number_of_trajectories : int
            The number of times to sample the chemical master equation. Each
            trajectory will be returned at the end of the simulation.
            Optional, defaults to 1.
        seed : int
            The random seed for the simulation. Optional, defaults to None.
        solver : gillespy.GillesPySolver
            The solver by which to simulate the model. This solver object may
            be initialized separately to specify an algorithm. Optional, 
            defaults to ssa solver.
        show_labels: bool (True)
            If true, simulation returns a list of trajectories, where each list entry is a dictionary containing key value pairs of species : trajectory.  If false, returns a numpy array with shape [traj_no, time, species]
        switch_tol: float
            Tolerance for Continuous/Stochastic representation of species, based on coefficient of variance for each step.
        tau_tol: float
            Relative error tolerance value for calculating tau step between 0.0 and 1.0
        integrator: String
            integrator to be used form scipy.integrate.ode. Options include 'vode', 'zvode', 'lsoda', 'dopri5', and 'dop835'.  For more details, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
        integrator_options: dictionary
            contains options to the scipy integrator. for a list of options, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html.
            Example use: {max_step : 0, rtol : .01}
        """
        if solver is not None:
            if ((isinstance(solver, type)
                    and issubclass(solver, GillesPySolver))) or issubclass(type(solver), GillesPySolver):
                solver_results = solver.run(model=self, t=self.tspan[-1], increment=self.tspan[-1] - self.tspan[-2], **solver_args)
            else:
                raise SimulationError(
                    "argument 'solver' to run() must be a subclass of GillesPySolver")
        else:
            from gillespy2.solvers.auto import SSASolver
            solver = SSASolver
            solver_results = SSASolver.run(model=self, t=self.tspan[-1],
                                      increment=self.tspan[-1] - self.tspan[-2], **solver_args)

        if isinstance(solver_results[0], (np.ndarray)):
            return solver_results

        if len(solver_results) is 1:
            return Results(data=solver_results[0], model=self, solver_name=solver.name)

        if len(solver_results) > 1:
            results_list = []
            for i in range(0,solver_args.get('number_of_trajectories')):
                results_list.append(Results(data=solver_results[i],model=self,solver_name=solver.name))
            return EnsembleResults(results_list)
        else:
            raise ValueError("number_of_trajectories must be non-negative and non-zero")


class Species(SortableObject):
    """
    Chemical species. Can be added to Model object to interact with other
    species or time.

    Attributes
    ----------
    name : str
        The name by which this species will be called in reactions and within
        the model.
    initial_value : int >= 0
        Initial population of this species. If this is not provided as an int,
        the type will be changed when it is added by numpy.int
    mode : str
        ***FOR USE WITH TauHybridSolver ONLY***
        Sets the mode of representation of this species for the TauHybridSolver,
        can be discrete, continuous, or dynamic.
        mode='dynamic' - Default, allows a species to be represented as
            either discrete or continuous
        mode='continuous' - Species will only be represented as continuous
        mode='discrete' - Species will only be represented as discrete
    """

    def __init__(self, name="", initial_value=0, mode='dynamic', allow_negative_populations=False):
        # A species has a name (string) and an initial value (positive integer)
        self.name = name
        self.mode = mode
        self.allow_negative_populations = allow_negative_populations

        mode_list = ['continuous', 'dynamic', 'discrete']
        if self.mode not in mode_list:
            raise SpeciesError('Species mode must be either \'continuous\', \'dynamic\', or \'discrete\'.')
        if mode == 'continuous':
            self.initial_value = np.float(initial_value)
        else:
            if not isinstance(initial_value, int): raise ValueError('Discrete values must be of type int.')
            self.initial_value = np.int(initial_value)
        if not allow_negative_populations:
            if self.initial_value < 0: raise ValueError('A species initial value must be \
non-negative unless allow_negative_populations=True')

    def __str__(self):
        return self.name


class Parameter(SortableObject):
    """
    A parameter can be given as an expression (function) or directly
    as a value (scalar). If given an expression, it should be
    understood as evaluable in the namespace of a parent Model.

    Attributes
    ----------
    name : str
        The name by which this parameter is called or referenced in reactions.
    expression : str
        String for a function calculating parameter values. Should be evaluable
        in namespace of Model.
    value : float
        Value of a parameter if it is not dependent on other Model entities.
    """

    def __init__(self, name="", expression=None, value=None):

        self.name = name
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.expression = expression
        if expression is not None:
            self.expression = str(expression)

        self.value = value

        # self.value is allowed to be None, but not self.expression. self.value
        # might not be evaluable in the namespace of this parameter, but defined
        # in the context of a model or reaction.
        if self.expression is None:
            raise TypeError

        if self.value is None:
            self.evaluate()

    def evaluate(self, namespace={}):
        """
        Evaluate the expression and return the (scalar) value in the given
        namespace.

        Attributes
        ----------
        namespace : dict (optional)
            The namespace in which to test evaluation of the parameter, if it
            involves other parameters, etc.
        """
        try:
            self.value = (float(eval(self.expression, namespace)))
        except:
            self.value = None

    def set_expression(self, expression):
        """
        Sets the expression for a parameter.
        """
        self.expression = expression
        # We allow expression to be passed in as a non-string type. Invalid
        # strings will be caught below. It is perfectly fine to give a scalar
        # value as the expression. This can then be evaluated in an empty
        # namespace to the scalar value.
        if expression is not None:
            self.expression = str(expression)

        if self.expression is None:
            raise TypeError

        self.evaluate()


class RateRule:
    def __init__(self, species=None, expression='', name=None):
        self.expression = expression
        self.species = species
        self.name = name



class Reaction(SortableObject):
    """
    Models a single reaction. A reaction has its own dicts of species
    (reactants and products) and parameters. The reaction's propensity
    function needs to be evaluable (and result in a non-negative scalar
    value) in the namespace defined by the union of those dicts.

    Attributes
    ----------
    name : str
        The name by which the reaction is called.
    reactants : dict
        The reactants that are consumed in the reaction, with stoichiometry. An
        example would be {R1 : 1, R2 : 2} if the reaction consumes two of R1 and
        one of R2, where R1 and R2 are Species objects.
    products : dict
        The species that are created by the reaction event, with stoichiometry.
        Same format as reactants.
    propensity_function : str
        The custom propensity fcn for the reaction. Must be evaluable in the
        namespace of the reaction using C operations.
    massaction : bool
        The switch to use a mass-action reaction. If set to True, a rate value
        is required.
    rate : float
        The rate of the mass-action reaction. Take care to note the units...
    annotation : str
        An optional note about the reaction.

    Notes
    ----------
    For a species that is NOT consumed in the reaction but is part of a mass
    action reaction, add it as both a reactant and a product.

    Mass-action reactions must also have a rate term added. Note that the input
    rate represents the mass-action constant rate independent of volume.
    """

    def __init__(self, name="", reactants={}, products={},
                 propensity_function=None, massaction=False,
                 rate=None, annotation=None):
        """
        Initializes the reaction using short-hand notation.
        """

        # Metadata
        self.name = name
        self.annotation = ""

        # We might use this flag in the future to automatically generate
        # the propensity function if set to True.
        if propensity_function is not None:
            self.massaction = False
            self.marate = None
        else:
            self.massaction = True

        self.propensity_function = propensity_function
        self.ode_propensity_function = propensity_function

        if self.propensity_function is not None and self.massaction:
            errmsg = ("Reaction {} You cannot set the propensity type to mass-action and simultaneously set a "
                      "propensity function.").format(self.name)
            raise ReactionError(errmsg)

        self.reactants = {}
        for r in reactants:
            rtype = type(r).__name__
            if rtype == 'instance':
                self.reactants[r.name] = reactants[r]
            else:
                self.reactants[r] = reactants[r]

        self.products = {}
        for p in products:
            rtype = type(p).__name__
            if rtype == 'instance':
                self.products[p.name] = products[p]
            else:
                self.products[p] = products[p]

        if self.massaction:
            self.type = "mass-action"
            if rate is None:
                self.marate = None
            else:
                self.marate = rate
                self.create_mass_action()
        else:
            self.type = "customized"

    def verify(self):
        """ Check if the reaction is properly formatted.
        Does nothing on sucesss, raises and error on failure."""
        if self.marate is None and self.propensity_function is None:
             raise ReactionError("You must specify either a mass-action rate or a propensity function")
        if len(self.reactants) == 0 and len(self.products) == 0:
            raise ReactionError("You must have a non-zero number of reactants or products.")

    def create_mass_action(self):
        """
        Initializes the mass action propensity function given
        self.reactants and a single parameter value.
        """
        # We support zeroth, first and second order propensities only.
        # There is no theoretical justification for higher order propensities.
        # Users can still create such propensities if they really want to,
        # but should then use a custom propensity.
        total_stoch = 0
        for r in sorted(self.reactants):
            total_stoch += self.reactants[r]
        if total_stoch > 2:
            raise ReactionError("Reaction: A mass-action reaction cannot involve more than two of one species or one "
                                "of two species.")
        # Case EmptySet -> Y

        propensity_function = self.marate.name
        ode_propensity_function = self.marate.name

        # There are only three ways to get 'total_stoch==2':
        for r in sorted(self.reactants):
            # Case 1: 2X -> Y
            if self.reactants[r] == 2:
                propensity_function = (propensity_function +
                                       "*" + str(r) + "*(" + str(r) + "-1)/vol")
                ode_propensity_function += '*' + str(r) + '*' + str(r)
            else:
                # Case 3: X1, X2 -> Y;
                propensity_function += "*" + str(r)
                ode_propensity_function += '*' + str(r)

        # Set the volume dependency based on order.
        order = len(self.reactants)
        if order == 2:
            propensity_function += "/vol"
        elif order == 0:
            propensity_function += "*vol"

        self.propensity_function = propensity_function
        self.ode_propensity_function = ode_propensity_function

    def setType(self, rxntype):
        """
        Sets reaction type to either "mass-action" or "customized"

        Attributes
        ----------
        rxntype : str
            Either "mass-action" or "customized"
        """
        if rxntype.lower() not in {'mass-action', 'customized'}:
            raise ReactionError("Invalid reaction type.")
        self.type = rxntype.lower()

        self.massaction = False if self.type == 'customized' else True

    def addReactant(self, S, stoichiometry):
        """
        Adds a reactant to the reaction (species that is consumed)

        Attributes
        ----------
        S : gillespy.Species
            Reactant to add to this reaction.
        stoichiometry : int
            The stoichiometry of the given reactant.
        """
        if stoichiometry <= 0:
            raise ReactionError("Reaction Stoichiometry must be a \
                                    positive integer.")
        self.reactants[S.name] = stoichiometry

    def addProduct(self, S, stoichiometry):
        """
        Adds a product to the reaction (species that is created)

        Attributes
        ----------
        S : gillespy.Species
            Product to add to this reaction.
        stoichiometry : int
            The stoichiometry of the given product.
        """
        self.products[S.name] = stoichiometry

    def Annotate(self, annotation):
        """
        Adds a note to the reaction

        Attributes
        ----------
        annotation : str
            An optional note about the reaction.
        """
        self.annotation = annotation

    def sanitized_propensity_function(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()))
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_propensity = self.propensity_function
        for id, name in enumerate(names):
            sanitized_propensity = sanitized_propensity.replace(name, "{"+str(id)+"}")
        return sanitized_propensity.format(*replacements)


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
            spec.append(md.species_to_element(model.listOfSpecies[sname]))
        md.document.append(spec)

        # Parameters
        params = eTree.Element('ParametersList')
        for pname in model.listOfParameters:
            params.append(md.parameter_to_element(
                model.listOfParameters[pname]))

        params.append(md.parameter_to_element(Parameter(name='vol', expression=model.volume)))

        md.document.append(params)

        # Reactions
        reacs = eTree.Element('ReactionsList')
        for rname in model.listOfReactions:
            reacs.append(md.reaction_to_element(model.listOfReactions[rname], model.volume))
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
        if model.name is "":
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
                print(model.volume)
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

    def species_to_element(self, S):
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

    def parameter_to_element(self, P):
        e = eTree.Element('Parameter')
        idElement = eTree.Element('Id')
        idElement.text = P.name
        e.append(idElement)
        expressionElement = eTree.Element('Expression')
        expressionElement.text = str(P.value)
        e.append(expressionElement)
        return e

    def reaction_to_element(self, R, model_volume):
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
