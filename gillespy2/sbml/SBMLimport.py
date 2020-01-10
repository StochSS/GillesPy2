import os
import gillespy2
import numpy
try:
    import libsbml
except ImportError:
    raise ImportError('libsbml is required to convert SBML files for GillesPy.')


init_state = {}

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

def __get_species(sbml_model, gillespy_model, errors):

    for i in range(sbml_model.getNumSpecies()):
        species = sbml_model.getSpecies(i)
        if species.getId() == 'EmptySet':
            errors.append([
                              "EmptySet species detected in model on line {0}. EmptySet is not an explicit species in "
                              "gillespy".format(
                                  species.getLine()), 0])
            continue
        name = species.getId()
        if species.isSetInitialAmount():
            int_value = int(species.getInitialAmount())
            value = species.getInitialAmount()
            if value == int_value:
                value = int_value
                mode = 'dynamic'
            else: mode = 'continuous'
        elif species.isSetInitialConcentration():
            value = species.getInitialConcentration()
            mode = 'continuous'
        else: # Assignment
            mode = 'continuous'
            value = 0

        constant = species.getConstant()
        boundary_condition = species.getBoundaryCondition()
        is_negative = value < 0.0
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

        gillespy_parameter = gillespy2.Parameter(name=name, expression=value)
        gillespy_model.add_parameter([gillespy_parameter])

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

def __get_reactions(sbml_model, gillespy_model):
    # local parameters
    for i in range(sbml_model.getNumReactions()):
        reaction = sbml_model.getReaction(i)
        kinetic_law = reaction.getKineticLaw()

        for j in range(kinetic_law.getNumParameters()):
            parameter = kinetic_law.getParameter(j)
            name = parameter.getId()
            value = parameter.getValue()
            gillespy_parameter = gillespy2.Parameter(name=name, expression=value)
            gillespy_model.add_parameter([gillespy_parameter])

    # reactions
    for i in range(sbml_model.getNumReactions()):
        reaction = sbml_model.getReaction(i)
        name = reaction.getId()

        reactants = {}
        products = {}

        for j in range(reaction.getNumReactants()):
            species = reaction.getReactant(j)

            if species.getSpecies() == "EmptySet":
                errors.append([
                                  "EmptySet species detected as reactant in reaction '{0}' on line {1}. EmptySet is "
                                  "not an explicit species in gillespy".format(
                                      reaction.getId(), species.getLine()), 0])
            else:
                reactants[species.getSpecies()] = species.getStoichiometry()

        # get products
        for j in range(reaction.getNumProducts()):
            species = reaction.getProduct(j)

            if species.getSpecies() == "EmptySet":
                errors.append([
                                  "EmptySet species detected as product in reaction '{0}' on line {1}. EmptySet is "
                                  "not an explicit species in gillespy".format(
                                      reaction.getId(), species.getLine()), 0])
            else:
                products[species.getSpecies()] = species.getStoichiometry()

        # propensity
        kinetic_law = reaction.getKineticLaw()
        propensity = kinetic_law.getFormula()

        gillespy_reaction = gillespy2.Reaction(name=name, reactants=reactants, products=products,
                                             propensity_function=propensity)

        gillespy_model.add_reaction([gillespy_reaction])

def __get_rules(sbml_model, gillespy_model, errors):
    for i in range(sbml_model.getNumRules()):
        rule = sbml_model.getRule(i)

        t = []
        
        '''
        if rule.isCompartmentVolume():
            t.append('compartment')
        if rule.isParameter():
            t.append('parameter')
        '''

        if rule.isAssignment():
            rule_name = rule.getId()
            rule_string = libsbml.formulaToL3String(rule.getMath())
            gillespy_rule = gillespy2.AssignmentRule(variable=rule_name,
                formula=rule_string)
            if rule_name in gillespy_model.listOfParameters:
                # Treat Non-Constant Parameters as Species
                value = gillespy_model.listOfParameters[rule_name].expression
                species = gillespy2.Species(name=rule_name,
                                            initial_value=value,
                                            mode='continuous')
                gillespy_model.delete_parameter(rule_name)
                gillespy_model.add_species([species])

            gillespy_model.add_assignment_rule(gillespy_rule)

        if rule.isRate():
            rule_name = rule.getId()
            rule_string = libsbml.formulaToL3String(rule.getMath())
            if rule_name in gillespy_model.listOfParameters:
                # Treat Non-Constant Parameters as Species
                value = gillespy_model.listOfParameters[rule_name].expression
                species = gillespy2.Species(name=rule_name,
                                            initial_value=value,
                                            mode='continuous')
                gillespy_model.delete_parameter(rule_name)
                gillespy_model.add_species([species])
           
            gillespy_rule = gillespy2.RateRule(species=gillespy_model.listOfSpecies[rule_name],
                expression=rule_string)
            gillespy_model.add_rate_rule(gillespy_rule)

        if rule.isAlgebraic():
            t.append('algebraic')

        if len(t) > 0:
            t[0] = t[0].capitalize()

            msg = ", ".join(t)
            msg += " rule"
        else:
            msg = "Rule"

        errors.append(["{0} '{1}' found on line '{2}' with equation '{3}'. gillespy does not support SBML Algebraic, compartment, or parameter Rules".format(
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
        function_string = libsbml.formulaToL3String(function.getMath())
        function_elements = function_string.replace('lambda(','')[:-1].split(', ')
        function_args = function_elements[:-1]
        function_function = function_elements[-1].replace('^', '**')
        gillespy_function = gillespy2.FunctionDefinition(name=function_name, function=function_function, args=function_args)
        gillespy_model.add_function_definition(gillespy_function)

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
            gillespy_assignment = gillespy2.EventAssignment(a.getVariable(),
                libsbml.formulaToL3String(a.getMath()))
            gillespy_assignments.append(gillespy_assignment)
        gillespy_event = gillespy2.Event(
            name=event.name, trigger=gillespy_trigger,
            assignments=gillespy_assignments, delay=delay,
            use_values_from_trigger_time=use_values_from_trigger_time)
        gillespy_model.add_event(gillespy_event)

def __get_initial_assignments(sbml_model, gillespy_model):
    for i in range(sbml_model.getNumInitialAssignments()):
        ia = sbml_model.getInitialAssignment(i)
        variable = ia.getId()
        expression = libsbml.formulaToL3String(ia.getMath())
        if variable in gillespy_model.listOfSpecies:
            print(variable)
            print(gillespy_model.listOfSpecies)
            print(init_state)
            gillespy_model.listOfSpecies[variable].initial_value = eval(expression, init_state)


def convert(filename, model_name=None, gillespy_model=None):

    sbml_model, errors = __read_sbml_model(filename)
    if model_name is None:
        model_name = sbml_model.getName()
    if gillespy_model is None:
        gillespy_model = gillespy2.Model(name=model_name)
    gillespy_model.units = "concentration"

    __get_species(sbml_model, gillespy_model, errors)
    __get_parameters(sbml_model, gillespy_model)
    __get_compartments(sbml_model, gillespy_model)
    __get_reactions(sbml_model, gillespy_model)
    __get_rules(sbml_model, gillespy_model, errors)
    __get_constraints(sbml_model, gillespy_model)
    __get_function_definitions(sbml_model, gillespy_model)
    __get_events(sbml_model, gillespy_model)
    __get_initial_assignments(sbml_model, gillespy_model)

    return gillespy_model, errors



