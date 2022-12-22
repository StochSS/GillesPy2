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

import uuid

from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.jsonify import Jsonify

class RateRule(SortableObject, Jsonify):
    """
    A RateRule is used to express equations that determine the rates of change
    of variables. This would correspond to a function in the form of dx/dt=f(W)

    :param name: Name of Rule
    :type name: str

    :param variable: Target Species/Parameter to be modified by rule
    :type variable: str

    :param formula: String representation of formula to be evaluated
    :type formula: str
    """

    def __init__(self, variable=None, formula='', name=None):
        if name in (None, ""):
            self.name = f'rr{uuid.uuid4()}'.replace('-', '_')
        else:
            self.name = name
        self.formula = formula
        self.variable = variable

    def __str__(self):
        var_name = self.variable if isinstance(self.variable, str) else self.variable.name
        return f"{self.name}: Var: {var_name}: {self.formula}"

    def sanitized_formula(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_formula = self.formula
        for i, name in enumerate(names):
            sanitized_formula = sanitized_formula.replace(name, "{" + str(i) + "}")
        return sanitized_formula.format(*replacements)
