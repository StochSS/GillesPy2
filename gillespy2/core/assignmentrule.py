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

from typing import (
    Dict,
    Union
)

from gillespy2.core.species import Species
from gillespy2.core.parameter import Parameter

from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.jsonify import Jsonify


class AssignmentRule(SortableObject, Jsonify):
    """
    An :class:`AssignmentRule` is used to express equations that set the values of
    variables. This would correspond to a function in the form of :code:`x = f(V)`.

    :param variable: Target :class:`~gillespy2.core.species.Species` or :class:`~gillespy2.core.parameter.Parameter`
        that this rule will modify.

    :param name: The name that this rule will be referenced by.
    :param formula: The string representation of the formula to evaluate.
    """

    def __init__(
        self, 
        variable: Union[Species, Parameter] = None, 
        formula: str = None, 
        name: str = None
    ):
        if name in (None, ""):
            self.name = f'ar{uuid.uuid4()}'.replace('-', '_')
        else:
            self.name = name
        self.variable = variable
        self.formula = formula

    def __str__(self) -> str:
        return self.variable + ': ' + self.formula

    def sanitized_formula(
        self, 
        species_mappings: Dict[str, str],
        parameter_mappings: Dict[str, str]
    ) -> str:
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_formula = self.formula
        for id, name in enumerate(names):
            sanitized_formula = sanitized_formula.replace(name, "{" + str(id) + "}")
        return sanitized_formula.format(*replacements)
