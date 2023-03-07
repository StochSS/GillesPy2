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
import unittest

from gillespy2.core.parameter import Parameter
from gillespy2.core.gillespyError import ParameterError

class TestParameter(unittest.TestCase):
    '''
    ################################################################################################
    Unit tests for gillespy2.Parameter.
    ################################################################################################
    '''
    def setUp(self):
        self.parameter = Parameter(name="test_parameter", expression="0.5")

    def test_constructor(self):
        """ Test the Parameter constructor. """
        parameter = Parameter(name="test_parameter", expression="0.5")
        self.assertEqual(parameter.name, "test_parameter")
        self.assertEqual(parameter.expression, "0.5")

    def test_constructor__no_name(self):
        """ Test the Parameter constructor without name. """
        with self.assertRaises(ParameterError):
            Parameter(expression="0.5")

    def test_constructor__invalid_name(self):
        """ Test the Parameter constructor with non-str name. """
        test_names = ["", None, 0, 0.5, [0]]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(ParameterError):
                    Parameter(name=test_name, expression="0.5")

    def test_constructor__no_expression(self):
        """ Test the Parameter constructor without expression. """
        with self.assertRaises(ParameterError):
            Parameter(name="test_parameter")

    def test_constructor__int_expression(self):
        """ Test the Parameter constructor with int expression. """
        parameter = Parameter(name="test_parameter", expression=1)
        self.assertEqual(parameter.expression, "1")

    def test_constructor__float_expression(self):
        """ Test the Parameter constructor with float expression. """
        parameter = Parameter(name="test_parameter", expression=0.5)
        self.assertEqual(parameter.expression, "0.5")

    def test_constructor__invaild_expression(self):
        """ Test the Parameter constructor with an invalid expression. """
        test_exps = [None, "", [2]]
        for test_exp in test_exps:
            with self.subTest(expression=test_exp):
                with self.assertRaises(ParameterError):
                    Parameter(name="test_name", expression=test_exp)

    def test___str__(self):
        """ Test Parameter.__str__ method. """
        self.assertIsInstance(str(self.parameter), str)

    def test__evaluate(self):
        """ Test Parameter._evaluate method. """
        self.parameter._evaluate()
        self.assertEqual(self.parameter.value, 0.5)

    def test__evaluate__int_expression(self):
        """ Test Parameter._evaluate method with int expression. """
        self.parameter.expression = 5
        self.parameter._evaluate()
        self.assertEqual(self.parameter.value, 5)

    def test__evaluate__float_expression(self):
        """ Test Parameter._evaluate method with float expression. """
        self.parameter.expression = 0.5
        self.parameter._evaluate()
        self.assertEqual(self.parameter.value, 0.5)

    def test__evaluate__improper_expression(self):
        """ Test Parameter._evaluate method with invalid expression. """
        self.parameter.expression = "[0.5]"
        with self.assertRaises(ParameterError):
            self.parameter._evaluate()

    def test__evaluate__invaild_expression(self):
        """ Test Parameter._evaluate with an invalid expression. """
        test_exps = [None, "", []]
        for test_exp in test_exps:
            with self.subTest(expression=test_exp):
                with self.assertRaises(ParameterError):
                    self.parameter.expression = test_exp
                    self.parameter._evaluate()

    def test__evaluate__parameter_in_namespace(self):
        """ Test Parameter._evaluate method with parameter in namespace. """
        self.parameter.expression = "k1 + 0.5"
        self.parameter._evaluate(namespace={"k1": 3})
        self.assertEqual(self.parameter.value, 3.5)

    def test__evaluate__species_in_namespace(self):
        """ Test Parameter._evaluate method with species in namespace. """
        self.parameter.expression = "S0 + 0.5"
        self.parameter._evaluate(namespace={"S0": 100})
        self.assertEqual(self.parameter.value, 100.5)

    def test__evaluate__component_not_in_namespace(self):
        """ Test Parameter._evaluate method with component missing from namespace. """
        test_comps = ["SO", "k1"]
        for test_comp in test_comps:
            with self.subTest(component=test_comp):
                self.parameter.expression = f"{test_comp} + 0.5"
                with self.assertRaises(ParameterError):
                    self.parameter._evaluate()

    def test_validate__invalid_name(self):
        """ Test Parameter.validate with non-str name. """
        test_names = ["", None, 0, 0.5, [0]]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(ParameterError):
                    self.parameter.name = test_name
                    self.parameter.validate()

    def test_validate__invaild_expression(self):
        """ Test Parameter.validate with an invalid expression. """
        test_exps = [None, "", [2]]
        for test_exp in test_exps:
            with self.subTest(expression=test_exp):
                with self.assertRaises(ParameterError):
                    self.parameter.expression = test_exp
                    self.parameter.validate()

    def test_comp_time_of_validate(self):
        """ Check the computation time of validate. """
        import time
        from datetime import datetime
        start = time.time()
        self.parameter.validate()
        tic = datetime.utcfromtimestamp(time.time() - start)
        print(f"Total time to run validate: {tic.strftime('%M mins %S secs %f msecs')}")
