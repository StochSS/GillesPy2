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
''' Unit tests module for gillespy2.AssignmentRule. '''
import re
import time
import unittest
from datetime import datetime

from gillespy2 import Species, Parameter, AssignmentRule
from gillespy2 import AssignmentRuleError

class TestAssignmentRule(unittest.TestCase):
    ''' Unit tests class for gillespy2.AssignmentRule. '''
    def setUp(self):
        self.assignmentrule = AssignmentRule(name="test_assignment_rule", variable="test_species", formula=29)

    def test_constructor(self):
        """ Test the AssignmentRule constructor. """
        assignmentrule = AssignmentRule(name="test_assignment_rule", variable="A", formula="sin(t)")
        self.assertEqual(assignmentrule.name, "test_assignment_rule")
        self.assertEqual(assignmentrule.variable, "A")
        self.assertEqual(assignmentrule.formula, "sin(t)")

    def test_constructor__no_name(self):
        """ Test the AssignmentRule constructor without name. """
        assignmentrule = AssignmentRule(variable="A", formula="sin(t)")
        self.assertIsNotNone(re.search("ar.*", assignmentrule.name))

    def test_constructor__name_is_none_or_empty(self):
        """ Test the AssignmentRule constructor with None or empty string as name. """
        test_names = [None, ""]
        for test_name in test_names:
            with self.subTest(name=test_name):
                assignmentrule = AssignmentRule(name=test_name, variable="A", formula="sin(t)")
                self.assertIsNotNone(re.search("ar.*", assignmentrule.name))

    def test_constructor__invalid_name(self):
        """ Test the AssignmentRule constructor with non-str name. """
        test_names = [0, 0.5, [0]]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(AssignmentRuleError):
                    AssignmentRule(name=test_name, variable="A", formula="sin(t)")

    def test_constructor__no_variable(self):
        """ Test the AssignmentRule constructor without variable. """
        with self.assertRaises(AssignmentRuleError):
            AssignmentRule(name="test_assignment_rule", formula="sin(t)")

    def test_constructor__species_variable(self):
        """ Test the AssignmentRule constructor with an species object as variable. """
        species = Species(name="s1", initial_value=5)
        assignmentrule = AssignmentRule(name="test_assignment_rule", variable=species, formula="sin(t)")
        self.assertIsInstance(assignmentrule.variable, str)
        self.assertEqual(assignmentrule.variable, "s1")

    def test_constructor__parameter_variable(self):
        """ Test the AssignmentRule constructor with an parameter object as variable. """
        parameter = Parameter(name="k1", expression="20")
        with self.assertRaises(AssignmentRuleError):
            AssignmentRule(name="test_assignment_rule", variable=parameter, formula="sin(t)")

    def test_constructor__invalid_variable(self):
        """ Test the AssignmentRule constructor with a variable that is an invalid type. """
        test_variables = [None, "", 5, 0.5, ["k1"]]
        for test_variable in test_variables:
            with self.subTest(variable=test_variable):
                with self.assertRaises(AssignmentRuleError):
                    AssignmentRule(name="test_assignment_rule", variable=test_variable, formula="sin(t)")

    def test_constructor__no_formula(self):
        """ Test the AssignmentRule constructor without formula. """
        with self.assertRaises(AssignmentRuleError):
            AssignmentRule(name="test_assignment_rule", variable="A")

    def test_constructor__int_formula(self):
        """ Test the AssignmentRule constructor with int formula. """
        assignmentrule = AssignmentRule(name="test_assignment_rule", variable="A", formula=5)
        self.assertEqual(assignmentrule.formula, "5")

    def test_constructor__float_formula(self):
        """ Test the AssignmentRule constructor with float formula. """
        assignmentrule = AssignmentRule(name="test_assignment_rule", variable="A", formula=0.5)
        self.assertEqual(assignmentrule.formula, "0.5")

    def test_constructor__invalid_formula(self):
        """ Test the AssignmentRule constructor with an invalid formula. """
        test_formulas = [None, "", [2]]
        for test_formula in test_formulas:
            with self.subTest(formula=test_formula):
                with self.assertRaises(AssignmentRuleError):
                    AssignmentRule(name="test_assignment_rule", variable="A", formula=test_formula)

    def test___str__(self):
        """ Test AssignmentRule.__str__ method. """
        self.assertIsInstance(str(self.assignmentrule), str)

    def test_validate__invalid_name(self):
        """ Test the AssignmentRule.validate with non-str name. """
        test_names = [None, "", 0, 0.5, [0]]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(AssignmentRuleError):
                    self.assignmentrule.name = test_name
                    self.assignmentrule.validate()

    def test_validate__parameter_variable(self):
        """ Test the AssignmentRule.validate with parameter variable. """
        parameter = Parameter(name="k1", expression="20")
        with self.assertRaises(AssignmentRuleError):
            self.assignmentrule.variable = parameter
            self.assignmentrule.validate()

    def test_validate__invalid_variable(self):
        """ Test the AssignmentRule.validate with an invalid variable. """
        test_variables = [None, "", 5, 0.5, ["k1"]]
        for test_variable in test_variables:
            with self.subTest(variable=test_variable):
                with self.assertRaises(AssignmentRuleError):
                    self.assignmentrule.variable = test_variable
                    self.assignmentrule.validate()

    def test_validate__invalid_formula(self):
        """ Test the AssignmentRule.validate with an invalid formula. """
        test_formulas = [None, "", 5, 0.5, [2]]
        for test_formula in test_formulas:
            with self.subTest(formula=test_formula):
                with self.assertRaises(AssignmentRuleError):
                    self.assignmentrule.formula = test_formula
                    self.assignmentrule.validate()

    def test_comp_time_of_validate(self):
        """ Check the computation time of validate. """
        start = time.time()
        self.assignmentrule.validate()
        tic = datetime.utcfromtimestamp(time.time() - start)
        msg = f"Total time to run validate on an assignment rule: {tic.strftime('%M mins %S secs %f msecs')}"
        print(f"\n<{'-'*88}>\n | {msg.ljust(84)} | \n<{'-'*88}>")
