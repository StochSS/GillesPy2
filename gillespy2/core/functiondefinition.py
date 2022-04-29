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
    List
)

from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.jsonify import Jsonify

class FunctionDefinition(SortableObject, Jsonify):
    """
    A :class:`FunctionDefinition` is the object representation of a mathematical function to be evaluated during the
    simulation of a :class:`~gillespy2.core.model.Model`.

    :param name: The name by which this function will be called and referenced.

    :param function: The mathematical function to be evaluated within the context of :code:`variables` and
        :class:`~gillespy2.core.model.Model` at simulation runtime.

    :param variables: A list of variable names to be used within the simulation context of :code:`function` at runtime.

    :raises TypeError: If :code:`function` is :code:`None`.
    """

    def __init__(
        self, 
        name: str = "", 
        function: str = None, 
        args: List[str] = []
    ):
        if function is None:
            raise TypeError("Function string provided for FunctionDefinition cannot be None")

        if name in (None, ""):
            self.name = f'fd{uuid.uuid4()}'.replace('-', '_')
        else:
            self.name = name
        self.function_string = function
        self.args = args


    def __str__(self) -> str:
        return f"{self.name}: Args: {self.args}, Expression: {self.function_string}"

    def get_arg_string(self) -> str:
        """
        Converts this :class:`FunctionDefinition` object's :code:`self.args` list into a comma-separated formatted
        string.

        :returns: `self.args` argument list as a comma-separated formatted string.
        """

        return ','.join(self.args)

    def sanitized_function(
        self, 
        species_mappings: Dict[str, str], 
        parameter_mappings: Dict[str, str]
    ) -> str:
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_function = self.function_string
        for id, name in enumerate(names):
            sanitized_function = sanitized_function.replace(name, "{" + str(id) + "}")
        return sanitized_function.format(*replacements)
