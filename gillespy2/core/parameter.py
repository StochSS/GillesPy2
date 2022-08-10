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

from gillespy2.core.jsonify import Jsonify
from gillespy2.core.sortableobject import SortableObject

from gillespy2.core.gillespyError import ParameterError

class Parameter(SortableObject, Jsonify):
    """
    A parameter can be given as an expression (function) or directly
    as a value (scalar). If given an expression, it should be
    understood as evaluable in the namespace of a parent Model.

    :param name: The name by which this parameter is called or referenced in \
                 reactions, rate rules, assignment rules, and events.
    :type name: str

    :param expression: String for a function calculating parameter values. Should be
                       evaluable in namespace of Model.
    :type expression: str

    :raises ParameterError: Arg is of invalid type.  Required arg set to None.  Arg value is outside of accepted bounds.
    """

    def __init__(self, name=None, expression=None):
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.value = None
        self.name = name

        if isinstance(expression, (int, float)):
            expression = str(expression)
        self.expression = expression

        self.validate()

    def __str__(self):
        return f"{self.name}: {self.expression}"

    def _evaluate(self, namespace=None):
        """
        Evaluate the expression and return the (scalar) value in the given
        namespace.

        :param namespace: The namespace in which to test evaulation of the parameter,
            if it involves other parameters, etc.
        :type namespace: dict

        :raises ParameterError: expression is of invalid type.  expression is set to None. \
                                expression is not evaluable within the given namespace.
        """
        if isinstance(self.expression, (int, float)):
            self.expression = str(self.expression)

        self.validate(coverage="expression")

        try:
            if namespace is None:
                namespace = {}
            self.value = float(eval(self.expression, namespace))
        except Exception as err:
            raise ParameterError(
                f"Could not evaluate expression: '{self.expression}'. Reason given: {err}."
            ) from err

    def sanitized_expression(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_expression = self.expression
        for i, name in enumerate(names):
            sanitized_expression = sanitized_expression.replace(
                name, "{"+str(i)+"}")
        return sanitized_expression.format(*replacements)

    def validate(self, expression=None, coverage="all"):
        """
        Validate the parameter.

        :param expression: String for a function calculating parameter values. Should be
                           evaluable in namespace of Model.
        :type expression: str

        :param coverage: The scope of attributes to validate.  Set to an attribute name to restrict validation \
                         to a specific attribute.
        :type coverage: str

        :raises ParameterError: Attribute is of invalid type.  Required attribute set to None.  \
                                Attribute value is outside of accepted bounds.
        """
        # Check name
        if coverage in ("all", "name"):
            if self.name is None:
                raise ParameterError("name can't be None type.")
            if not isinstance(self.name, str):
                raise ParameterError("name must be of type str.")
            if self.name == "":
                raise ParameterError("name can't be an empty string.")

        # Check expression
        if coverage in ("all", "expression"):
            if expression is None:
                expression = self.expression

            if expression is None:
                raise ParameterError("initial_value can't be None type.")
            if not isinstance(expression, str):
                raise ParameterError("expression must be of type str, float, or int.")
            if expression == "":
                raise ParameterError("expression can't be an empty string.")
