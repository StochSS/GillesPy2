import sys
import numpy as np
from .gillespy2 import example_models


def interior_trigonometric(token, model):
    token_list = token.split("(", 1)
    # We can assume that the first token is our trigonometric function.
    # I tend to alias Numpy to np as per convention, ask PI if cool.
    trig_dict = {"sin": "numpy.sine",
                 "cos": "numpy.cos",
                 "tan": "numpy.tan"}
    # I'll add more later, I don't see them crop up a lot in my parameter list.
    return "{}{}".format(trig_dict[token_list[0]], parenthetical_concatenation(token_list[1], model))


def interior_exponentiation(token, model):
    token_list = token.split("(", 1)
    return "(**{}".format(parenthetical_concatenation(token_list[1], model))


def interior_comparator(self, token):
    # Ask about this it could get weird.
    pass


def interior_conjunction(self, token):
    pass


def interior_piecewise(self, token):
    pass


def parenthetical_concatenation(token, model):
    return "({})".format(function_rules_parser(token[1:-1], model))


def function_rules_parser(string, model):
    operand_list = ['+', '-', '*', '/']
    # Working List of functions found as parameter.
    # Ask about conditional statements. I've seen these in matlab models, but they are present.
    function_dict = {'sin': interior_trigonometric,
                     'pow': interior_exponentiation,
                     'exp': interior_exponentiation,
                     'gt': interior_comparator,
                     'geq': interior_comparator,
                     'lt': interior_comparator,
                     'leq': interior_comparator,
                     'eq': interior_comparator,
                     'and': interior_conjunction,
                     'piecewise': interior_piecewise}
    # Approximately 8% of all BioDatabase Parameter rules are floats or integers that can simply be coerced.
    try:
        return float(string)
    except (TypeError, ValueError):
        # Causal observation shows that most of the string is white space delimited save the interior of parenthetical
        # Statements. Will have to have interior loop on tokens to further separate. Possible recursion?
        tokens = string.split(" ")
        ret_expression = ""
        for token in tokens:
            if token[0] == '(' or token[0] == ')':
                ret_expression+=token[0]
                token = token[1:]
                tokens.insert(0, token)
                continue
            if token in operand_list:
                ret_expression += token
                continue
            if token in function_dict.keys():
                function_designation = token[:2]
            #If any of these tokens are not in these lists we can assume they refer to a species of some type or parameter
            else:
                #Mixed case of numbers and variables
                #REFP. (Sorry.)
                try:
                    x = float(token)
                    ret_expression += token
                except (TypeError, ValueError):
                    try:
                        x = int(token)
                        ret_expression += token
                    except (TypeError, ValueError):
                        if token in model.listOfSpecies.keys():
                            ret_expression += "model.listOfSpecies['{}'].initial_value".format(token)
                        else:
                            print("Could not parse string token: {}, please check model.".format(token))
        return ret_expression


model = Trichloroethylene()
testString = "2.0 + (TCE + Epoxide)"
x =(function_rules_parser(testString, model))
print(x)
eval(x)
