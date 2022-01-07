"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.gillespyError import *
from gillespy2.core.jsonify import Jsonify

class Parameter(SortableObject, Jsonify):
    """
    A parameter can be given as an expression (function) or directly
    as a value (scalar). If given an expression, it should be
    understood as evaluable in the namespace of a parent Model.

    :param name: The name by which this parameter is called or referenced in reactions.
    :type name: str

    :param expression: String for a function calculating parameter values. Should be
        evaluable in namespace of Model.
    :type expression: str
    """

    def __init__(self, name="", expression=None):

        self.value = None
        self.name = name
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.

        if expression is None:
            raise ParameterError("Parameter expression can not be none")

        self.expression = str(expression)

    def __str__(self):
        return self.name + ': ' + str(self.expression)

    def set_expression(self, expression):
        """
        Sets the expression for a parameter.
        """
        # We allow expression to be passed in as a non-string type. Invalid
        # strings will be caught below. It is perfectly fine to give a scalar
        # value as the expression. This can then be evaluated in an empty
        # namespace to the scalar value.
        from gillespy2.core import log

        log.warning("'Parameter.set_expression' has been deprecated, future versions of GillesPy2 will not support"
                    " this function. To set expression within a parameter, use Parameter.expression = expression")

        if expression is None:
            raise ParameterError("Parameter expression can not be none")

        self.expression = str(expression)



    def _evaluate(self, namespace={}):
        """
        Evaluate the expression and return the (scalar) value in the given
        namespace.

        :param namespace: The namespace in which to test evaulation of the parameter,
            if it involves other parameters, etc.
        :type namespace: dict
        """

        try:
            self.value = (float(eval(self.expression, namespace)))
        except Exception as error:
            raise ParameterError("Could not evaluate expression: {}.".format(str(error))) from error

    def sanitized_expression(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_expression = self.expression
        for id, name in enumerate(names):
            sanitized_expression = sanitized_expression.replace(
                name, "{"+str(id)+"}")
        return sanitized_expression.format(*replacements)
