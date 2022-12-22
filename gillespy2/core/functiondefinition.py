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

from gillespy2.core.gillespyError import FunctionDefinitionError

class FunctionDefinition(SortableObject, Jsonify):
    """
    Object representation defining an evaluable function to be used during
    simulation of a GillesPy2 model

    :param name: Name of the function to be made and called
    :type name: str

    :param function: Defined function body of operation to be performed.
    :type function: str

    :param variables: String names of Variables to be used as arguments to function.
    :type variables: list[str]
    """

    def __init__(self, name="", function=None, args=[]):
        if function is None:
            raise FunctionDefinitionError("Function string provided for FunctionDefinition cannot be None")

        if name in (None, ""):
            self.name = f'fd{uuid.uuid4()}'.replace('-', '_')
        else:
            self.name = name
        self.function_string = function
        self.args = args


    def __str__(self):
        return f"{self.name}: Args: {self.args}, Expression: {self.function_string}"

    def get_arg_string(self) -> str:
        """
        Convert function's argument list into a comma-separated formatted string.

        :returns: Argument list as a comma-separated formatted string.
        """
        return ','.join(self.args)

    def sanitized_function(self, species_mappings, parameter_mappings):
        '''
        Sanitize the function definition function.

        :param species_mappings: Mapping of species names to sanitized species names.
        :type species_mappings: dict

        :param parameter_mappings: Mapping of parameter names to sanitized parameter names.
        :type parameter_mappings: dict

        :returns: The sanitized function.
        :rtype: str
        '''
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_function = self.function_string
        for i, name in enumerate(names):
            sanitized_function = sanitized_function.replace(name, "{" + str(i) + "}")
        return sanitized_function.format(*replacements)
