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
''' Unit tests module for gillespy2.RateRule. '''
import re
import time
import unittest
from datetime import datetime

from gillespy2 import Species, Parameter, RateRule
from gillespy2 import RateRuleError

class TestParameter(unittest.TestCase):
    ''' Unit tests class for gillespy2.RateRule. '''
    def setUp(self):
        self.raterule = RateRule(name="test_rate_rule", variable="test_species", formula=29)

    def test_constructor(self):
        """ Test the RateRule constructor. """
        raterule = RateRule(name="test_rate_rule", variable="A", formula="sin(t)")
        self.assertEqual(raterule.name, "test_rate_rule")
        self.assertEqual(raterule.variable, "A")
        self.assertEqual(raterule.formula, "sin(t)")

    def test_constructor__no_name(self):
        """ Test the RateRule constructor without name. """
        raterule = RateRule(variable="A", formula="sin(t)")
        self.assertIsNotNone(re.search("rr.*", raterule.name))

    def test_constructor__name_is_none_or_empty(self):
        """ Test the RateRule constructor with None or empty string as name. """
        test_names = [None, ""]
        for test_name in test_names:
            with self.subTest(name=test_name):
                raterule = RateRule(name=test_name, variable="A", formula="sin(t)")
                self.assertIsNotNone(re.search("rr.*", raterule.name))

    def test_constructor__invalid_name(self):
        """ Test the RateRule constructor with non-str name. """
        test_names = [0, 0.5, [0]]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(RateRuleError):
                    RateRule(name=test_name, variable="A", formula="sin(t)")

    def test_constructor__no_variable(self):
        """ Test the RateRule constructor without variable. """
        with self.assertRaises(RateRuleError):
            RateRule(name="test_rate_rule", formula="sin(t)")

    def test_constructor__species_variable(self):
        """ Test the RateRule constructor with an species object as variable. """
        species = Species(name="s1", initial_value=5)
        raterule = RateRule(name="test_rate_rule", variable=species, formula="sin(t)")
        self.assertIsInstance(raterule.variable, str)
        self.assertEqual(raterule.variable, "s1")

    def test_constructor__parameter_variable(self):
        """ Test the RateRule constructor with an parameter object as variable. """
        parameter = Parameter(name="k1", expression="20")
        with self.assertRaises(RateRuleError):
            RateRule(name="test_rate_rule", variable=parameter, formula="sin(t)")

    def test_constructor__invalid_variable(self):
        """ Test the RateRule constructor with a variable that is an invalid type. """
        test_variables = [None, "", 5, 0.5, ["k1"]]
        for test_variable in test_variables:
            with self.subTest(variable=test_variable):
                with self.assertRaises(RateRuleError):
                    RateRule(name="test_rate_rule", variable=test_variable, formula="sin(t)")

    def test_constructor__no_formula(self):
        """ Test the RateRule constructor without formula. """
        with self.assertRaises(RateRuleError):
            RateRule(name="test_rate_rule", variable="A")

    def test_constructor__int_formula(self):
        """ Test the RateRule constructor with int formula. """
        raterule = RateRule(name="test_rate_rule", variable="A", formula=5)
        self.assertEqual(raterule.formula, "5")

    def test_constructor__float_formula(self):
        """ Test the RateRule constructor with float formula. """
        raterule = RateRule(name="test_rate_rule", variable="A", formula=0.5)
        self.assertEqual(raterule.formula, "0.5")

    def test_constructor__invalid_formula(self):
        """ Test the RateRule constructor with an invalid formula. """
        test_formulas = [None, "", [2]]
        for test_formula in test_formulas:
            with self.subTest(formula=test_formula):
                with self.assertRaises(RateRuleError):
                    RateRule(name="test_rate_rule", variable="A", formula=test_formula)

    def test___str__(self):
        """ Test RateRule.__str__ method. """
        self.assertIsInstance(str(self.raterule), str)

    def test_validate__invalid_name(self):
        """ Test the RateRule.validate with non-str name. """
        test_names = [None, "", 0, 0.5, [0]]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(RateRuleError):
                    self.raterule.name = test_name
                    self.raterule.validate()

    def test_validate__parameter_variable(self):
        """ Test the RateRule.validate with parameter variable. """
        parameter = Parameter(name="k1", expression="20")
        with self.assertRaises(RateRuleError):
            self.raterule.variable = parameter
            self.raterule.validate()

    def test_validate__invalid_variable(self):
        """ Test the RateRule.validate with an invalid variable. """
        test_variables = [None, "", 5, 0.5, ["k1"]]
        for test_variable in test_variables:
            with self.subTest(variable=test_variable):
                with self.assertRaises(RateRuleError):
                    self.raterule.variable = test_variable
                    self.raterule.validate()

    def test_validate__invalid_formula(self):
        """ Test the RateRule.validate with an invalid formula. """
        test_formulas = [None, "", 5, 0.5, [2]]
        for test_formula in test_formulas:
            with self.subTest(formula=test_formula):
                with self.assertRaises(RateRuleError):
                    self.raterule.formula = test_formula
                    self.raterule.validate()

    def test_comp_time_of_validate(self):
        """ Check the computation time of validate. """
        start = time.time()
        self.raterule.validate()
        tic = datetime.utcfromtimestamp(time.time() - start)
        msg = f"Total time to run validate on a rate rule: {tic.strftime('%M mins %S secs %f msecs')}"
        print(f"\n<{'-'*88}>\n | {msg.ljust(84)} | \n<{'-'*88}>")
