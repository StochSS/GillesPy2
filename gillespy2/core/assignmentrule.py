# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2023 GillesPy2 developers.

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

from gillespy2.core.species import Species
from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.jsonify import Jsonify

from gillespy2.core.gillespyError import AssignmentRuleError

class AssignmentRule(SortableObject, Jsonify):
    """
    An AssignmentRule is used to express equations that set the values of
    variables.  This would correspond to a function in the form of x = f(V)

    :param name: Name of the Rule
    :type name: str

    :param variable: Target Species to be modified by rule
    :type variable: str

    :param formula: String representation of formula to be evaluated
    :type formula: str
    """

    def __init__(self, name=None, variable=None, formula=None):
        if name in (None, ""):
            name = f'ar{uuid.uuid4()}'.replace('-', '_')
        if isinstance(formula, (int, float)):
            formula = str(formula)

        self.name = name
        self.formula = formula

        self.validate(variable=variable)

        if variable is not None:
            vtype = type(variable).__name__
            if vtype == 'Species':
                variable = variable.name
        self.variable = variable

    def __str__(self):
        var_name = self.variable if isinstance(self.variable, str) else self.variable.name
        return f"{self.name}: Var: {var_name}: {self.formula}"

    def _create_sanitized_assignment_rule(self, n_ndx, species_mappings, parameter_mappings):
        name = f"AR{n_ndx}"
        variable = species_mappings[self.variable.name]
        formula = self.sanitized_formula(species_mappings, parameter_mappings)
        return AssignmentRule(name=name, formula=formula, variable=variable)

    def sanitized_formula(self, species_mappings, parameter_mappings):
        '''
        Sanitize the assignment rule formula.

        :param species_mappings: Mapping of species names to sanitized species names.
        :type species_mappings: dict

        :param parameter_mappings: Mapping of parameter names to sanitized parameter names.
        :type parameter_mappings: dict

        :returns: The sanitized formula.
        :rtype: str
        '''
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_formula = self.formula
        for i, name in enumerate(names):
            sanitized_formula = sanitized_formula.replace(name, "{" + str(i) + "}")
        return sanitized_formula.format(*replacements)

    def validate(self, variable=None, coverage="all"):
        """
        Validate the assignment rule.

        :param variable: Target Species to be modified by rule
        :type variable: str

        :param coverage: The scope of attributes to validate.  Set to an attribute name to restrict validation \
                         to a specific attribute.
        :type coverage: str

        :raises SpeciesError: Attribute is of invalid type.  Required attribute set to None.  \
                              Attribute is value outside of accepted bounds.
        """
        # Check name
        if coverage in ("all", "name"):
            if self.name is None:
                raise AssignmentRuleError("name can't be None type.")
            if not isinstance(self.name, str):
                raise AssignmentRuleError(f"name must be of type str not {type(self.name)}.")
            if self.name == "":
                raise AssignmentRuleError("name can't be an empty string.")

        # Check variable
        if coverage in ("all", "variable"):
            if variable is None:
                if not hasattr(self, "variable") or self.variable is None:
                    raise AssignmentRuleError("Assignment rules must have a variable.")
                variable = self.variable

            if not (isinstance(variable, (str, Species)) or type(variable).__name__ == 'Species'):
                raise AssignmentRuleError("variable must be of type str or GillesPy2.Species.")
            if variable == "":
                raise AssignmentRuleError("variable can't be an empty string.")

        # Check formula
        if coverage in ("all", "formula"):
            if not isinstance(self.formula, str):
                raise AssignmentRuleError("formula must be of type str.")
            if self.formula == "":
                raise AssignmentRuleError("formula can't be an empty string.")
