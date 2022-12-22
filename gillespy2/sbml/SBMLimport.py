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
import re
import math

try:
    import libsbml
except ImportError as err:
    raise ImportError('libsbml is required to convert SBML files for GillesPy.') from err

import numpy as np

import gillespy2
from gillespy2.core.gillespyError import InvalidModelError, SBMLError

init_state = {'INF': np.inf, 'NaN': np.nan}
postponed_evals = {}
eval_globals = math.__dict__.copy()
def piecewise(*args):
    args = list(args)
    sol = None
    if len(args) % 2:
        args.append(True)
    for i, arg in enumerate(args):
        if not i % 2:
            continue
        if arg:
            sol = args[i-1]
            break
    return sol
eval_globals['piecewise'] = piecewise

def __read_sbml_model(filename):

    document = libsbml.readSBML(filename)
    errors = []
    error_count = document.getNumErrors()
    if error_count > 0:
        for i in range(error_count):
            error = document.getError(i)
            converter_code = 0
            converter_code = -10

            errmsg = f"SBML {error.getSeverityAsString()}, code {error.getErrorId()}, "
            errmsg += f"line {error.getLine()}: {error.getMessage()}"
            errors.append([errmsg, converter_code])
    if min([code for error, code in errors] + [0]) < 0:
        return None, errors
    sbml_model = document.getModel()

    return sbml_model, errors

def __get_math(formula):
    math_str = libsbml.formulaToL3String(formula)
    replacements = {
        r'\bln\b': 'log',
        r'\^': '**',
        r'\&\&': 'and',
        r'\|\|': 'or'
        }
    for old, new in replacements.items():
        math_str = re.sub(old, new, math_str)
    return math_str

def __get_species(sbml_model, gillespy_model, errors):

    for i in range(sbml_model.getNumSpecies()):
        species = sbml_model.getSpecies(i)
        name = species.getId()
        population = None
        concentration = None
        value = None
        mode = 'continuous'
        cid = sbml_model.getCompartment(species.getCompartment()).getName()
        if species.isSetInitialAmount(): # Amount is provided
            population = species.getInitialAmount()
        if species.isSetInitialConcentration(): # Concentration is provided
            concentration = species.getInitialConcentration()
        if population is None and concentration is None: # Assignment
            value = 0
        elif species.getHasOnlySubstanceUnits(): # Treat as population
            if population is not None: # If population is provided
                value = population
                int_value = int(species.getInitialAmount())
                if population == int_value:
                    population = int_value
            else: # Else convert concentration to population
                postponed_evals[name] = f'{name} * {cid}'
                value = concentration
        else: # Treat as concentration
            if concentration is not None: # If concentration is provided
                value = concentration
            else: # Else convert population to concentration
                postponed_evals[name] = f'{name} / {cid}'
                value = population

        constant = species.getConstant()
        boundary_condition = species.getBoundaryCondition()
        is_negative = value < 0
        gillespy_species = gillespy2.Species(name=name, initial_value=value,
                                                allow_negative_populations=is_negative, mode=mode,
                                                constant=constant, boundary_condition=boundary_condition)
        gillespy_model.add_species([gillespy_species])
        init_state[name] = value

def __get_parameters(sbml_model, gillespy_model):

    for i in range(sbml_model.getNumParameters()):
        parameter = sbml_model.getParameter(i)
        name = parameter.getId()
        if parameter.isSetValue():
            value = parameter.getValue()
        else:
            value = 0
        init_state[name] = value

        # GillesPy2 represents non-constant parameters as species
        if parameter.isSetConstant():
            gillespy_parameter = gillespy2.Parameter(name=name, expression=value)
            gillespy_model.add_parameter([gillespy_parameter])
        else:
            gillespy_species = gillespy2.Species(name=name,initial_value=value)
            gillespy_model.add_species([gillespy_species])

def __get_compartments(sbml_model, gillespy_model):
    for i in range(sbml_model.getNumCompartments()):
        compartment = sbml_model.getCompartment(i)
        name = compartment.getId()
        if compartment.isSetSize():
            value = compartment.getSize()
        else:
            value = 1

        if name == "vol":
            gillespy_model.volume = value
        else:
            gillespy_parameter = gillespy2.Parameter(name=name, expression=value)
            init_state[name] = value
            gillespy_model.add_parameter([gillespy_parameter])

def __traverse_math(node, old_id, new_id):
    if node is None:
        return
    for i in range(node.getNumChildren()):
        if node.getChild(i).getName() == old_id:
            new_node = libsbml.ASTNode()
            new_node.setName(new_id)
            node.replaceChild(i, new_node)
        __traverse_math(node.getChild(i), old_id, new_id)

def __get_kinetic_law(sbml_model, gillespy_model, reaction):
    kinetic_law = reaction.getKineticLaw()

    if kinetic_law is None:
        raise InvalidModelError(
            f"Failed to load SBML model: Reaction '{reaction}' is missing its propensity function."
        )

    tree = kinetic_law.getMath()
    params = kinetic_law.getListOfParameters()
    local_params = kinetic_law.getListOfLocalParameters()
    for i in range(kinetic_law.getNumLocalParameters()):
        local_param = local_params.get(i)
        old_id = local_param.getId()
        new_id = (f'{reaction.getId()}_{local_param.getId()}')
        __traverse_math(tree, old_id, new_id)
        local_param.setId(new_id)
        gillespy_parameter = gillespy2.Parameter(name=new_id, expression=local_param.getValue())
        gillespy_model.add_parameter([gillespy_parameter])
    for i in range(kinetic_law.getNumParameters()):
        param = params.get(i)
        if not param.getId() in gillespy_model.listOfParameters:
            gillespy_parameter = gillespy2.Parameter(name=param.getId(), expression=param.getValue())
            gillespy_model.add_parameter([gillespy_parameter])
    return tree

def __get_reactions(sbml_model, gillespy_model, errors):
    # reactions
    for i in range(sbml_model.getNumReactions()):
        reaction = sbml_model.getReaction(i)
        name = reaction.getId()
        tree = __get_kinetic_law(sbml_model, gillespy_model, reaction)
        propensity = __get_math(tree)
        reactants = {}
        products = {}

        r_set = set()
        p_set = set()
        # get reactants
        for j in range(reaction.getNumReactants()):
            reactant = reaction.getReactant(j)
            species = reactant.getSpecies()

            if species == "EmptySet":
                continue
            stoichiometry = reactant.getStoichiometry()
            if isinstance(stoichiometry, float):
                if int(stoichiometry) != stoichiometry:
                    logmsg = f"Reaction {name} contains a float stoichiometry for reactant {species}.  "
                    logmsg += "Please check your model as this may cause inaccuracies in the results."
                    gillespy2.log.warning(logmsg)
                stoichiometry = int(stoichiometry)

            if species in r_set:
                reactants[species] += stoichiometry
            else:
                r_set.add(species)
                reactants[species] = stoichiometry

        # get products
        for j in range(reaction.getNumProducts()):
            product = reaction.getProduct(j)
            species = product.getSpecies()

            if species == "EmptySet":
                continue
            stoichiometry = product.getStoichiometry()
            if isinstance(stoichiometry, float):
                if int(stoichiometry) != stoichiometry:
                    logmsg = f"Reaction {name} contains a float stoichiometry for product {species}.  "
                    logmsg += "Please check your model as this may cause inaccuracies in the results."
                    gillespy2.log.warning(logmsg)
                stoichiometry = int(stoichiometry)

            if species in p_set:
                products[species] += stoichiometry
            else:
                p_set.add(species)
                products[species] = stoichiometry

        gillespy_reaction = gillespy2.Reaction(name=name, reactants=reactants, products=products,
                                             propensity_function=propensity)

        gillespy_model.add_reaction([gillespy_reaction])

def __get_rules(sbml_model, gillespy_model, errors):
    for i in range(sbml_model.getNumRules()):
        rule = sbml_model.getRule(i)

        # If the SBML object does not contain an ID attribute create a unique rule_name from the variable name.
        rule_name = rule.getIdAttribute()
        if rule_name == "":
            rule_name = f"ar_{rule.getId()}"

        rule_variable = rule.getVariable()
        rule_string = __get_math(rule.getMath())
        if rule_variable in gillespy_model.listOfParameters:
            # Treat Non-Constant Parameters as Species
            value = gillespy_model.listOfParameters[rule_variable].expression
            species = gillespy2.Species(name=rule_variable,
                                        initial_value=value,
                                        allow_negative_populations=True,
                                        mode='continuous')
            gillespy_model.delete_parameter(rule_variable)
            gillespy_model.add_species([species])

        t = []

        if rule.isAssignment():
            _ = eval(rule_string, {**eval_globals, **init_state})
            postponed_evals[rule_variable] = rule_string
            gillespy_rule = gillespy2.AssignmentRule(name=rule_name, variable=rule_variable,
                formula=rule_string)
            gillespy_model.add_assignment_rule(gillespy_rule)
            init_state[gillespy_rule.variable.name]=eval(gillespy_rule.formula, {**init_state, **eval_globals})

        if rule.isRate():
            gillespy_rule = gillespy2.RateRule(name=rule_name, variable=rule_variable,
                formula=rule_string)
            gillespy_model.add_rate_rule(gillespy_rule)

        if rule.isAlgebraic():
            t.append('algebraic')

            if len(t) > 0:
                t[0] = t[0].capitalize()

                msg = ", ".join(t)
                msg += " rule"
            else:
                msg = "Rule"

            errmsg = f"{msg} '{rule.getId()}' found on line '{rule.getLine()}' with equation "
            errmsg += f"'{libsbml.formulaToString(rule.getMath())}'. gillespy does not support SBML Algebraic Rules."
            errors.append([errmsg, -5])

def __get_constraints(sbml_model, gillespy_model, errors):
    for i in range(sbml_model.getNumConstraints()):
        constraint = sbml_model.getConstraint(i)

        errmsg = f"Constraint '{constraint.getId()}' found on line '{constraint.getLine()}' with equation "
        errmsg += f"'{libsbml.formulaToString(constraint.getMath())}'. gillespy does not support SBML Constraints"
        errors.append([errmsg, -5])

def __get_function_definitions(sbml_model, gillespy_model):
    # TODO:
    # DOES NOT CURRENTLY SUPPORT ALL MATHML
    # ALSO DOES NOT SUPPORT NON-MATHML
    for i in range(sbml_model.getNumFunctionDefinitions()):
        function = sbml_model.getFunctionDefinition(i)
        function_name = function.getId()
        function_tree = function.getMath()
        num_nodes = function_tree.getNumChildren()
        function_args = [function_tree.getChild(i).getName() for i in range(num_nodes-1)]
        function_string = __get_math(function_tree.getChild(num_nodes-1))
        gillespy_function = gillespy2.FunctionDefinition(
            name=function_name, function=function_string, args=function_args
        )
        gillespy_model.add_function_definition(gillespy_function)
        init_state[gillespy_function.name] = gillespy_function.function_string

def __get_events(sbml_model, gillespy_model):
    for i in range(sbml_model.getNumEvents()):
        event = sbml_model.getEvent(i)
        gillespy_assignments = []

        trigger = event.getTrigger()
        delay = event.getDelay()
        if delay is not None:
            delay = libsbml.formulaToL3String(delay.getMath())
        expression=libsbml.formulaToL3String(trigger.getMath())
        expression = expression.replace('&&', ' and ').replace('||', ' or ')
        initial_value = trigger.getInitialValue()
        persistent = trigger.getPersistent()
        use_values_from_trigger_time = event.getUseValuesFromTriggerTime()
        gillespy_trigger = gillespy2.EventTrigger(
            expression=expression, initial_value=initial_value, persistent=persistent
        )
        assignments = event.getListOfEventAssignments()
        for assign in assignments:
            # Convert Non-Constant Parameter to Species
            if assign.getVariable() in gillespy_model.listOfParameters:
                gillespy_species = gillespy2.Species(
                    name=assign.getVariable(),
                    initial_value=gillespy_model.listOfParameters[assign.getVariable()].expression,
                    mode='continuous', allow_negative_populations=True
                )
                gillespy_model.delete_parameter(assign.getVariable())
                gillespy_model.add_species([gillespy_species])

            gillespy_assignment = gillespy2.EventAssignment(variable=assign.getVariable(),
                expression=__get_math(assign.getMath()))
            gillespy_assignments.append(gillespy_assignment)
        gillespy_event = gillespy2.Event(
            name=event.getId(), trigger=gillespy_trigger,
            assignments=gillespy_assignments, delay=delay,
            use_values_from_trigger_time=use_values_from_trigger_time)
        gillespy_model.add_event(gillespy_event)

def __get_initial_assignments(sbml_model, gillespy_model):

    for i in range(sbml_model.getNumInitialAssignments()):
        init_assign = sbml_model.getInitialAssignment(i)
        variable = init_assign.getId()
        expression = __get_math(init_assign.getMath())
        assigned_value = eval(expression, {**init_state, **eval_globals})
        init_state[variable] = assigned_value
        if assigned_value != assigned_value:
            assigned_value = expression
            postponed_evals[variable] = expression

        if variable in gillespy_model.listOfSpecies:
            gillespy_model.listOfSpecies[variable].initial_value = assigned_value
        elif variable in gillespy_model.listOfParameters:
            gillespy_model.listOfParameters[variable].expression = assigned_value

def __resolve_evals(gillespy_model, init_state):
    while True:
        successful = []
        if len(postponed_evals):
            for var, expr in postponed_evals.items():
                try:
                    assigned_value = eval(expr, {**eval_globals, **init_state})
                except Exception:
                    assigned_value = np.nan
                if assigned_value == assigned_value:
                    successful.append(var)
                    init_state[var] = assigned_value
                    if var in gillespy_model.listOfSpecies:
                        gillespy_model.listOfSpecies[var].initial_value = assigned_value
                    elif var in gillespy_model.listOfParameters:
                        gillespy_model.listOfParameters[var].value = assigned_value
        if len(successful) <= 0:
            break
        for var in successful:
            del postponed_evals[var]

def convert(filename, model_name=None, gillespy_model=None, report_silently_with_sbml_error=False):
    sbml_model, errors = __read_sbml_model(filename)

    if sbml_model is None:
        if report_silently_with_sbml_error:
            return None, errors
        errs = '\n\t'.join(errors)
        raise SBMLError(f"SBML model import failed.  Reason Given: \n\t{errs}")

    if len(errors) > 0 and not report_silently_with_sbml_error:
        from gillespy2 import log # pylint: disable=import-outside-toplevel
        errs = '\n\t'.join(errors)
        errmsg = f"Error were detected in the SBML model.  Error: \n\t{errs}"
        log.warning(errmsg)

    if model_name is None:
        model_name = sbml_model.getName()
    if gillespy_model is None:
        gillespy_model = gillespy2.Model(name=model_name)
    gillespy_model.units = "concentration"

    __get_function_definitions(sbml_model, gillespy_model)
    __get_parameters(sbml_model, gillespy_model)
    __get_species(sbml_model, gillespy_model, errors)
    __get_compartments(sbml_model, gillespy_model)
    __get_reactions(sbml_model, gillespy_model, errors)
    __get_rules(sbml_model, gillespy_model, errors)
    __get_constraints(sbml_model, gillespy_model, errors)
    __get_events(sbml_model, gillespy_model)
    __get_initial_assignments(sbml_model, gillespy_model)
    __resolve_evals(gillespy_model, init_state)

    return gillespy_model, errors
