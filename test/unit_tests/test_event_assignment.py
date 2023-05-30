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
import time
import unittest
from datetime import datetime

from gillespy2 import Species, Parameter, EventAssignment
from gillespy2 import EventAssignmentError

class TestEventAssignment(unittest.TestCase):
    ''' Unit tests class for gillespy2.EventAssignment. '''
    def setUp(self):
        self.eventassignment = EventAssignment(variable="event_flag", expression="1")

    def test_constructor(self):
        """ Test the EventAssignment constructor. """
        eventassignment = EventAssignment(variable="event_flag", expression="1")
        self.assertEqual(eventassignment.variable, "event_flag")
        self.assertEqual(eventassignment.expression, "1")

    def test_constructor__no_variable(self):
        """ Test the EventAssignment constructor without variable. """
        with self.assertRaises(EventAssignmentError):
            EventAssignment(expression="1")

    def test_constructor__species_variable(self):
        """ Test the EventAssignment constructor with an species object as variable. """
        species = Species(name="s1", initial_value=5)
        eventassignment = EventAssignment(variable=species, expression="1")
        self.assertIsInstance(eventassignment.variable, str)
        self.assertEqual(eventassignment.variable, "s1")

    def test_constructor__parameter_variable(self):
        """ Test the EventAssignment constructor with an parameter object as variable. """
        parameter = Parameter(name="k1", expression="20")
        with self.assertRaises(EventAssignmentError):
            EventAssignment(variable=parameter, expression="1")

    def test_constructor__invalid_variable(self):
        """ Test the EventAssignment constructor with a variable that is an invalid type. """
        test_variables = [None, "", 5, 0.5, ["k1"]]
        for test_variable in test_variables:
            with self.subTest(variable=test_variable):
                with self.assertRaises(EventAssignmentError):
                    EventAssignment(variable=test_variable, expression="1")

    def test_constructor__no_expression(self):
        """ Test the EventAssignment constructor without expression. """
        with self.assertRaises(EventAssignmentError):
            EventAssignment(variable="event_flag")

    def test_constructor__int_expression(self):
        """ Test the EventAssignment constructor with int expression. """
        eventassignment = EventAssignment(variable="event_flag", expression=1)
        self.assertEqual(eventassignment.expression, "1")

    def test_constructor__float_expression(self):
        """ Test the EventAssignment constructor with float expression. """
        eventassignment = EventAssignment(variable="event_flag", expression=0.5)
        self.assertEqual(eventassignment.expression, "0.5")

    def test_constructor__invalid_expression(self):
        """ Test the EventAssignment constructor with an invalid expression. """
        test_expressions = [None, "", [2]]
        for test_expression in test_expressions:
            with self.subTest(expression=test_expression):
                with self.assertRaises(EventAssignmentError):
                    EventAssignment(variable="event_flag", expression=test_expression)

    def test___str__(self):
        """ Test EventAssignment.__str__ method. """
        self.assertIsInstance(str(self.eventassignment), str)

    def test_validate__parameter_variable(self):
        """ Test the EventAssignment.validate with parameter variable. """
        parameter = Parameter(name="k1", expression="20")
        with self.assertRaises(EventAssignmentError):
            self.eventassignment.variable = parameter
            self.eventassignment.validate()

    def test_validate__invalid_variable(self):
        """ Test the EventAssignment.validate with an invalid variable. """
        test_variables = [None, "", 5, 0.5, ["k1"]]
        for test_variable in test_variables:
            with self.subTest(variable=test_variable):
                with self.assertRaises(EventAssignmentError):
                    self.eventassignment.variable = test_variable
                    self.eventassignment.validate()

    def test_validate__invalid_formula(self):
        """ Test the EventAssignment.validate with an invalid expression. """
        test_expressions = [None, "", 5, 0.5, [2]]
        for test_expression in test_expressions:
            with self.subTest(expression=test_expression):
                with self.assertRaises(EventAssignmentError):
                    self.eventassignment.expression = test_expression
                    self.eventassignment.validate()

    def test_comp_time_of_validate(self):
        """ Check the computation time of validate. """
        start = time.time()
        self.eventassignment.validate()
        tic = datetime.utcfromtimestamp(time.time() - start)
        msg = f"Total time to run validate on an event assignment: {tic.strftime('%M mins %S secs %f msecs')}"
        print(f"\n<{'-'*88}>\n | {msg.ljust(84)} | \n<{'-'*88}>")
