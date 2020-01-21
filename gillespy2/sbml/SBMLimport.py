'''
Convert a model from an SBML file into a GillesPy2 model.
'''

import os
import gillespy2
from   gillespy2 import Model, Parameter, Species, Reaction
import numpy

#
# Internal constants.
#

_UNSUPPORTED_FEATURE = -5
_CONVERSION_ERROR = -10

#
# Exported functions.
#

def convert(filename, model_name=None, model=None):
    try:
        import libsbml
    except ImportError:
        raise ImportError('Failed to import python-libsbml -- cannot convert SBML model')

    # FIXME test if file exists
    document = libsbml.readSBML(filename)

    # Return early if the SBML model contains errors.
    errors = []
    if document.getNumErrors() > 0:
        for i in range(document.getNumErrors()):
            error = document.getError(i)
            errors.append(["SBML {0}, error #{1}, line {2}: {3}".format(
                error.getSeverityAsString(), error.getErrorId(),
                error.getLine(), error.getMessage()), _CONVERSION_ERROR])
    if any(code < 0 for error, code in errors):
        return None, errors

    # The SBML model is free of basic errors.  Extract the components of the
    # model into a more convenient form for our work later below.

    sbml_model     = document.getModel()
    comp_list      = [sbml_model.getCompartment(i) for i in range(sbml_model.getNumCompartments())]
    species_list   = [sbml_model.getSpecies(i) for i in range(sbml_model.getNumSpecies())]
    parameter_list = [sbml_model.getParameter(i) for i in range(sbml_model.getNumParameters())]
    reaction_list  = [sbml_model.getReaction(i) for i in range(sbml_model.getNumReactions())]
    event_list     = [sbml_model.getEvent(i) for i in range(sbml_model.getNumEvents())]

    # We have a bit of a chicken-and-egg issue: we need to create a model
    # object in order to add add components like species, but we need to know
    # the volume first, and we can't get the volume without reading the
    # compartments in which the species are located.  Additional wrinkle: we
    # currently don't handle models with more than 1 reaction volume, but
    # it's possible for a model to define more than one compartment and yet
    # only use one as the volume for all species.  So, we iterate over all the
    # species first and check if more than one is actually used.

    if len(comp_list) > 1:
        used_compartment_names = set(s.getCompartment() for s in species)
        if len(used_compartment_names) > 1:
            msg = ('{} compartments found in model.'.format(len(comp_list))
                   + ' GillesPy2 only supports a single reaction volume.')
            errors.append([msg, _UNSUPPORTED_FEATURE])
            return None, errors
        else:
            # Multiple compartments, but all species are placed in only one.
            used_compartment = sbml_model.getCompartment(used_compartment_names[0])
            volume = used_compartment.getVolume()
    else:
        # There's only one compartment.
        volume = comp_list[0].getVolume()

    # If we made it this far, start building up the GillesPy2 model.

    if model_name is None:
        model_name = sbml_model.getName()
    if model is None:
        model = Model(name=model_name, volume=volume)

    # FIXME
    model.units = "concentration"

    for species in species_list:
        name = species.getId()
        if species.isSetInitialAmount():
            int_value = int(species.getInitialAmount())
            value = species.getInitialAmount()
            if value == int_value:
                value = int_value
                mode = 'dynamic'
            else:
                mode = 'continuous'
        elif species.isSetInitialConcentration():
            value = species.getInitialConcentration()
            mode = 'continuous'
        else:
            rule = sbml_model.getRule(species.getId())
            if rule:
                msg = ""
                if rule.isAssignment():
                    msg = "assignment "
                elif rule.isRate():
                    msg = "rate "
                elif rule.isAlgebraic():
                    msg = "algebraic "

                msg += "rule"

                errors.append([
                                  "Species '{0}' does not have any initial conditions. Associated {1} '{2}' found, "
                                  "but {1}s are not supported in gillespy. Assuming initial condition 0".format(
                                      species.getId(), msg, rule.getId()), 0])
            else:
                errors.append([
                                  "Species '{0}' does not have any initial conditions or rules. Assuming initial "
                                  "condition 0".format(
                                      species.getId()), 0])

            value = 0

        is_negative = value < 0.0
        is_boundary = species.getBoundaryCondition()
        is_constant = species.getConstant()

        model.add_species(Species(name=name, initial_value=value, mode=mode,
                                  allow_negative_populations=is_negative,
                                  boundary=is_boundary, constant=is_constant))

    # Global parameters.
    for parameter in parameter_list:
        name = parameter.getId()
        value = parameter.getValue()
        model.add_parameter(Parameter(name=name, expression=value))

    # Local parameters.
    # FIXME name collisions.
    for reaction in reaction_list:
        kl = reaction.getKineticLaw()
        for parameter in [kl.getParameter(i) for i in range(kl.getNumParameters())]:
            name = parameter.getId()
            value = parameter.getValue()
            model.add_parameter(Parameter(name=name, expression=value))

    # We treat compartments like parameters.
    for compartment in comp_list:
        name = compartment.getId()
        value = compartment.getSize()
        model.add_parameter(Parameter(name=name, expression=value))

    # Reactions.
    for reaction in reaction_list:
        reactants = {}
        products = {}
        for i in range(reaction.getNumReactants()):
            species_ref = reaction.getReactant(i)
            reactants[species_ref.getSpecies()] = species_ref.getStoichiometry()
        for i in range(reaction.getNumProducts()):
            species_ref = reaction.getProduct(i)
            products[species_ref.getSpecies()] = species_ref.getStoichiometry()

        propensity = reaction.getKineticLaw().getFormula()
        model.add_reaction(Reaction(name=reaction.getId(), reactants=reactants,
                                    products=products, propensity_function=propensity))

    # Events.
    for event in event_list:
        gillespy_assignments = []
        trigger = event.getTrigger()
        delay = event.getDelay()
        if delay is not None:
            delay = libsbml.formulaToL3String(delay.getMath())
        expression = libsbml.formulaToL3String(trigger.getMath())
        expression = expression.replace('&&', ' and ').replace('||', ' or ')
        initial_value = trigger.getInitialValue()
        persistent = trigger.getPersistent()
        use_values_from_trigger_time = event.getUseValuesFromTriggerTime()
        gillespy_trigger = gillespy2.EventTrigger(expression=expression,
            initial_value=initial_value, persistent=persistent)
        for a in event.getListOfEventAssignments():
            gillespy_assignment = gillespy2.EventAssignment(a.getVariable(),
                libsbml.formulaToL3String(a.getMath()))
            gillespy_assignments.append(gillespy_assignment)
        gillespy_event = gillespy2.Event(
            name=event.name, trigger=gillespy_trigger,
            assignments=gillespy_assignments, delay=delay,
            use_values_from_trigger_time=use_values_from_trigger_time)
        model.add_event(gillespy_event)

    # Report things that are currently unsupported.
    for i in range(sbml_model.getNumRules()):
        rule = sbml_model.getRule(i)

        t = []

        if rule.isCompartmentVolume():
            t.append('compartment')
        if rule.isParameter():
            t.append('parameter')
        elif rule.isAssignment():
            t.append('assignment')
        elif rule.isRate():
            t.append('rate')
        elif rule.isAlgebraic():
            t.append('algebraic')

        if len(t) > 0:
            t[0] = t[0].capitalize()

            msg = ", ".join(t)
            msg += " rule"
        else:
            msg = "Rule"

        errors.append(["{0} '{1}' found on line '{2}' with equation '{3}'. gillespy does not support SBML Rules".format(
            msg, rule.getId(), rule.getLine(), libsbml.formulaToString(rule.getMath())), -5])

    for i in range(sbml_model.getNumConstraints()):
        constraint = sbml_model.getConstraint(i)
        errors.append([
                          "Constraint '{0}' found on line '{1}' with equation '{2}'. gillespy does not support SBML "
                          "Constraints".format(
                              constraint.getId(), constraint.getLine(), libsbml.formulaToString(constraint.getMath())),
                          -5])

    for i in range(sbml_model.getNumFunctionDefinitions()):
        function = sbml_model.getFunctionDefinition(i)
        errors.append([
                          "Function '{0}' found on line '{1}' with equation '{2}'. gillespy does not support SBML "
                          "Function Definitions".format(
                              function.getId(), function.getLine(), libsbml.formulaToString(function.getMath())), -5])

    return model, errors
