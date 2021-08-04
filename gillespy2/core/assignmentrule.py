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
from gillespy2.core.jsonify import Jsonify


class AssignmentRule(SortableObject, Jsonify):
    """
    An AssignmentRule is used to express equations that set the values of
    variables.  This would correspond to a function in the form of x = f(V)

    :param name: Name of the Rule
    :type name: str

    :param variable: Target Species/Parameter to be modified by rule
    :type variable: str

    :param formula: String representation of formula to be evaluated
    :type formula: str
    """

    def __init__(self, variable=None, formula=None, name=None):
        self.variable = variable
        self.formula = formula
        self.name = name

    def __str__(self):
        return self.variable + ': ' + self.formula

    def sanitized_formula(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_formula = self.formula
        for id, name in enumerate(names):
            sanitized_formula = sanitized_formula.replace(name, "{" + str(id) + "}")
        return sanitized_formula.format(*replacements)