import os
import gillespy2
import numpy as np
import math
try:
    import libsbml
except ImportError:
    raise ImportError('libsbml is required to convert SBML files for GillesPy.')


init_state = {'INF': np.inf, 'NaN': np.nan}
postponed_evals = {}
eval_globals = math.__dict__.copy()
def piecewise(*args):
    args = list(args)
    sol = None
    if len(args) % 2: args.append(True)
    for i, arg in enumerate(args):
        if not i % 2: continue
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

            errors.append(["SBML {0}, code {1}, line {2}: {3}".format(error.getSeverityAsString(), error.getErrorId(),
                                                                      error.getLine(), error.getMessage()),
                           converter_code])
    if min([code for error, code in errors] + [0]) < 0:
        return None, errors
    sbml_model = document.getModel()

    return sbml_model, errors

def __get_math(math):
    math_str = libsbml.formulaToL3String(math)
    replacements = {
        'ln': 'log',
        '^': '**',
        '&&': 'and',
        '||': 'or'
        }
    for old, new in replacements.items():
        math_str = math_str.replace(old, new)
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
                postponed_evals[name] = '{} * {}'.format(name, cid)
                value = concentration
        else: # Treat as concentration
            if concentration is not None: # If concentration is provided
                value = concentration
            else: # Else convert population to concentration
                postponed_evals[name] = '{} / {}'.format(name, cid)
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
        value = parameter.getValue()
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
        value = compartment.getSize()

        gillespy_parameter = gillespy2.Parameter(name=name, expression=value)
        init_state[name] = value
        gillespy_model.add_parameter([gillespy_parameter])

    '''
    for i in range(sbml_model.getNumCompartments()):
        compartment = sbml_model.getCompartment(i)
        vol = compartment.getSize()
        gillespy_model.volume = vol

        errors.append([
                          "Compartment '{0}' found on line '{1}' with volume '{2}' and dimension '{3}'. gillespy "
                          "assumes a single well-mixed, reaction volume".format(
                              compartment.getId(), compartment.getLine(), compartment.getVolume(),
                              compartment.getSpatialDimensions()), -5])
    '''
def traverse_math(node, old_id, new_id):
    if node is None: return
    for i in range(node.getNumChildren()):
        if node.getChild(i).getName() == old_id:
            new_node = libsbml.ASTNode()
            new_node.setName(new_id)
            node.replaceChild(i, new_node)
        traverse_math(node.getChild(i), old_id, new_id)

def __get_kinetic_law(sbml_model, gillespy_model, reaction):
    kinetic_law = reaction.getKineticLaw()
    tree = kinetic_law.getMath()
    params = kinetic_law.getListOfParameters()
    local_params = kinetic_law.getListOfLocalParameters()
    for i in range(kinetic_law.getNumLocalParameters()):
        lp = local_params.get(i)
        old_id = lp.getId()
        new_id = ('{}_{}'.format(reaction.getId(), lp.getId()))
        traverse_math(tree, old_id, new_id)
        lp.setId(new_id)
        gillespy_parameter = gillespy2.Parameter(name=new_id, expression=lp.getValue())
        gillespy_model.add_parameter([gillespy_parameter])
    for i in range(kinetic_law.getNumParameters()):
        p = params.get(i)
        if not p.getId() in gillespy_model.listOfParameters:
            gillespy_parameter = gillespy2.Parameter(name=p.getId(), expression=p.getValue())
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
            species = reaction.getReactant(j)
            if species.getSpecies() == "EmptySet": continue
            else:
                if species.getSpecies() in r_set:
                    reactants[species.getSpecies()] += species.getStoichiometry()
                else:
                    r_set.add(species.getSpecies())
                    reactants[species.getSpecies()] = species.getStoichiometry()

        # get products
        for j in range(reaction.getNumProducts()):
            species = reaction.getProduct(j)

            if species.getSpecies() == "EmptySet": continue
            else:
                if species.getSpecies() in p_set:
                    products[species.getSpecies()] += species.getStoichiometry()
                else:
                    p_set.add(species.getSpecies())
                    products[species.getSpecies()] = species.getStoichiometry()

        gillespy_reaction = gillespy2.Reaction(name=name, reactants=reactants, products=products,
                                             propensity_function=propensity)

        gillespy_model.add_reaction([gillespy_reaction])

def __get_rules(sbml_model, gillespy_model, errors):
    for i in range(sbml_model.getNumRules()):
        rule = sbml_model.getRule(i)
        rule_name = rule.getId()
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
            assign_value = eval(rule_string, {**eval_globals, **init_state})
            postponed_evals[rule_variable] = rule_string
            gillespy_rule = gillespy2.AssignmentRule(name=rule_name, variable=rule_variable,
                formula=rule_string)
            gillespy_model.add_assignment_rule(gillespy_rule)
            init_state[gillespy_rule.variable]=eval(gillespy_rule.formula, {**init_state, **eval_globals})

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

            errors.append(["{0} '{1}' found on line '{2}' with equation '{3}'. gillespy does not support SBML Algebraic Rules".format(
                msg, rule.getId(), rule.getLine(), libsbml.formulaToString(rule.getMath())), -5])

def __get_constraints(sbml_model, gillespy_model):
    for i in range(sbml_model.getNumConstraints()):
        constraint = sbml_model.getConstraint(i)

        errors.append([
                          "Constraint '{0}' found on line '{1}' with equation '{2}'. gillespy does not support SBML "
                          "Constraints".format(
                              constraint.getId(), constraint.getLine(), libsbml.formulaToString(constraint.getMath())),
                          -5])

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
        gillespy_function = gillespy2.FunctionDefinition(name=function_name, function=function_string, args=function_args)
        gillespy_model.add_function_definition(gillespy_function)
        init_state[gillespy_function.name] = gillespy_function.function

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
        gillespy_trigger = gillespy2.EventTrigger(expression=expression, 
            initial_value=initial_value, persistent=persistent)
        assignments = event.getListOfEventAssignments()
        for a in assignments:
            # Convert Non-Constant Parameter to Species
            if a.getVariable() in gillespy_model.listOfParameters:
                gillespy_species = gillespy2.Species(name=a.getVariable(),
                                                        initial_value=gillespy_model.listOfParameters[a.getVariable()].expression,
                                                        mode='continuous', allow_negative_populations=True)
                gillespy_model.delete_parameter(a.getVariable())
                gillespy_model.add_species([gillespy_species])

            gillespy_assignment = gillespy2.EventAssignment(a.getVariable(),
                __get_math(a.getMath()))
            gillespy_assignments.append(gillespy_assignment)
        gillespy_event = gillespy2.Event(
            name=event.getId(), trigger=gillespy_trigger,
            assignments=gillespy_assignments, delay=delay,
            use_values_from_trigger_time=use_values_from_trigger_time)
        gillespy_model.add_event(gillespy_event)

def __get_initial_assignments(sbml_model, gillespy_model):

    for i in range(sbml_model.getNumInitialAssignments()):
        ia = sbml_model.getInitialAssignment(i)
        variable = ia.getId()
        expression = __get_math(ia.getMath())
        assigned_value = eval(expression, {**init_state, **eval_globals})
        init_state[variable] = assigned_value
        if assigned_value != assigned_value:
            assigned_value = expression
            postponed_evals[variable] = expression

        if variable in gillespy_model.listOfSpecies:
            gillespy_model.listOfSpecies[variable].initial_value = assigned_value
        elif variable in gillespy_model.listOfParameters:
            gillespy_model.listOfParameters[variable].set_expression(assigned_value)

def __resolve_evals(gillespy_model, init_state):
    while True:
        successful = []
        if len(postponed_evals):
            for var, expr in postponed_evals.items():
                try: assigned_value = eval(expr, {**eval_globals, **init_state})
                except: assigned_value = np.nan
                if assigned_value == assigned_value:
                    successful.append(var)
                    init_state[var] = assigned_value
                    if var in gillespy_model.listOfSpecies:
                        gillespy_model.listOfSpecies[var].initial_value = assigned_value
                    elif var in gillespy_model.listOfParameters:
                        gillespy_model.listOfParameters[var].value = assigned_value
        if not len(successful): break
        for var in successful: del postponed_evals[var]

def convert(filename, model_name=None, gillespy_model=None):

    sbml_model, errors = __read_sbml_model(filename)
    if sbml_model is None:
        return None, errors
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
    __get_constraints(sbml_model, gillespy_model)
    __get_events(sbml_model, gillespy_model)
    __get_initial_assignments(sbml_model, gillespy_model)
    __resolve_evals(gillespy_model, init_state)


    return gillespy_model, errors



