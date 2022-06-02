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

import ast
import json

from gillespy2.core.jsonify import ComplexJsonCoder

def __add_events(model, events):
    for name, event in events.items():
        s_event = {"compID":model['defaultID'],
                   "name": name,
                   "annotation": "",
                   "delay": event.delay,
                   "priority": event.priority,
                   "triggerExpression": event.trigger.expression,
                   "initialValue": event.trigger.value,
                   "persistent": event.trigger.persistent,
                   "useValuesFromTriggerTime": event.use_values_from_trigger_time,
                   "eventAssignments": []}

        __add_event_assignments(model=model, event=s_event, assignments=event.assignments)

        model['eventsCollection'].append(s_event)
        model['defaultID'] += 1


def __add_event_assignments(model, event, assignments):
    for assignment in assignments:
        name = assignment.variable.name
        try:
            variable = __get_species(species=model['species'], name=name)
        except IndexError:
            variable = __get_parameter(parameters=model['parameters'], name=name)

        s_assignment = {"variable": variable,
                        "expression": assignment.expression}
        event['eventAssignments'].append(s_assignment)


def __add_function_definitions(model, function_definitions):
    for name, function_definition in function_definitions.items():
        variables = function_definition.get_arg_string()
        expression = function_definition.function_string
        function = "lambda({0}, {1})".format(variables, expression)
        signature = "{0}({1})".format(name, variables)

        s_function_definition = {"compID":model['defaultID'],
                                 "name":name,
                                 "function":function,
                                 "expression":expression,
                                 "variables":variables,
                                 "signature":signature,
                                 "annotation": ""}
        model['functionDefinitions'].append(s_function_definition)
        model['defaultID'] += 1


def __add_parameters(model, parameters):
    for name, parameter in parameters.items():
        try:
            expression = ast.literal_eval(parameter.expression)
        except:
            expression = parameter.expression
        s_parameter = {"compID":model['defaultID'],
                       "name":name,
                       "expression":str(expression),
                       "annotation": ""}
        model['parameters'].append(s_parameter)
        model['defaultID'] += 1


def __add_reactions(model, reactions):
    for name, reaction in reactions.items():
        s_reaction = {"compID":model['defaultID'],
                      "name":name,
                      "reactionType": "custom-propensity",
                      "massaction": False,
                      "propensity": reaction.propensity_function,
                      "annotation": "",
                      "rate": {},
                      "types": [],
                      "reactants": [],
                      "products": []}

        for key in ['reactants', 'products']:
            __add_stoich_species(s_reaction=s_reaction, reaction=reaction,
                                                key=key, species=model['species'])
        __add_summary(reaction=s_reaction)

        model['reactions'].append(s_reaction)
        model['defaultID'] += 1


def __add_rules(model, r_type, rules):
    for name, rule in rules.items():
        try:
            variable = __get_species(species=model['species'], name=rule.variable.name)
        except IndexError:
            variable = __get_parameter(parameters=model['parameters'], name=rule.variable.name)

        s_rule = {"compID":model['defaultID'],
                  "name":name,
                  "expression":rule.formula,
                  "type":r_type,
                  "variable":variable,
                  "annotation": ""}
        model['rules'].append(s_rule)
        model['defaultID'] += 1


def __add_species(model, species):
    modes = []

    for name, specie in species.items():
        if specie.mode is not None and specie.mode not in modes:
            modes.append(specie.mode)
        s_species = {"compID":model['defaultID'],
                     "name":name,
                     "value":specie.initial_value,
                     "mode":specie.mode,
                     "switchTol": specie.switch_tol,
                     "switchMin": specie.switch_min,
                     "isSwitchTol": specie.switch_min == 0,
                     "annotation": "",
                     "diffusionConst":0,
                     "types": []}
        model['species'].append(s_species)
        model['defaultID'] += 1

    if not modes:
        model['defaultMode'] = ""
    elif len(modes) > 1:
        model['defaultMode'] = "dynamic"
    else:
        model['defaultMode'] = modes[0]


def __add_stoich_species(s_reaction, reaction, key, species):
    source = reaction.reactants if key == "reactants" else reaction.products
    for specie, ratio in source.items():
        stoich_species = {"ratio":ratio,
                          "specie":__get_species(species=species, name=specie.name)}
        s_reaction[key].append(stoich_species)


def __add_summary(reaction):
    r_summary = __build_element(reaction['reactants'])
    p_summary = __build_element(reaction['products'])
    reaction['summary'] = f"{r_summary} \\rightarrow {p_summary}"


def __build_element(stoich_species):
    if not stoich_species:
        return "\\emptyset"

    elements = []
    for species in stoich_species:
        name = species['specie']['name']
        ratio = species['ratio']
        element = f"{ratio}{name}" if ratio > 1 else name
        elements.append(element)
    return '+'.join(elements)


def __get_parameter(parameters, name):
    return list(filter(lambda parameter: parameter['name'] == name, parameters))[0]


def __get_species(species, name):
    return list(filter(lambda specie: specie['name'] == name, species))[0]


def __write_to_file(model, path):
    with open(path, "w") as model_file:
        json.dump(model, model_file, indent=4, sort_keys=True, cls=ComplexJsonCoder)


def export(model, path=None, return_stochss_model=False):
    """
    GillesPy model to StochSS converter

    :param model: GillesPy model to be converted to StochSS
    :type model: gillespy.Model

    :param path: Path to the StochSS file for conversion
    :type path: str
    """

    model.compile_prep()

    if path is None:
        path = f"{model.name}.mdl"

    if model.tspan is None:
        model_settings = {"endSim": 20, "timeStep": 0.05}
    else:
        model_settings = {
            "endSim": model.tspan[-1],
            "timeStep": model.tspan[-1] / len(model.tspan)
        }

    s_model = {"is_spatial": False,
               "defaultID": 1,
               "annotation": "",
               "volume": model.volume,
               "modelSettings": model_settings,
                "species": [],
                "initialConditions": [],
                "parameters": [],
                "reactions": [],
                "rules": [],
                "eventsCollection": [],
                "functionDefinitions": []
              }

    __add_species(model=s_model, species=model.get_all_species())
    __add_parameters(model=s_model, parameters=model.get_all_parameters())
    __add_reactions(model=s_model, reactions=model.get_all_reactions())
    __add_events(model=s_model, events=model.get_all_events())
    __add_rules(model=s_model, r_type='Rate Rule', rules=model.get_all_rate_rules())
    __add_rules(model=s_model, r_type='Assignment Rule', rules=model.get_all_assignment_rules())
    __add_function_definitions(model=s_model,
                               function_definitions=model.get_all_function_definitions())

    if return_stochss_model:
        return s_model

    __write_to_file(model=s_model, path=path)
    return path
