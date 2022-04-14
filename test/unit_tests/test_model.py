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
import sys
sys.path.insert(1, "../")
import unittest

from example_models import RobustModel
from gillespy2.core.gillespyError import *

class TestModel(unittest.TestCase):

    def setUp(self):
        self.model = RobustModel()

    def test_model_add__assignment_rule(self):
        from gillespy2 import AssignmentRule
        ar1 = AssignmentRule(name="ar1", variable="k1", formula="29")
        self.model.add(ar1)
        self.assertIn("ar1", self.model.listOfAssignmentRules)

    def test_model_add__event(self):
        from gillespy2 import Event, EventTrigger, EventAssignment
        ea = EventAssignment(name="ea", variable="k1", expression="29")
        et = EventTrigger(expression="t > 29")
        e2 = Event(name="e2", trigger=et, assignments=[ea])
        self.model.add(e2)
        self.assertIn("e2", self.model.listOfEvents)

    def test_model_add__function_definition(self):
        from gillespy2 import FunctionDefinition
        divide = FunctionDefinition(name="divide", function="x / y", args=["x", "y"])
        self.model.add(divide)
        self.assertIn("divide", self.model.listOfFunctionDefinitions)

    def test_model_add__invalid_component(self):
        with self.assertRaises(ModelError):
            self.model.add("29")

    def test_model_add__parameter(self):
        from gillespy2 import Parameter
        k3 = Parameter(name="k3", expression=29)
        self.model.add(k3)
        self.assertIn("k3", self.model.listOfParameters)

    def test_model_add__rate_rule(self):
        from gillespy2 import RateRule
        rr3 = RateRule(name="rr3", variable="k1", formula="29")
        self.model.add(rr3)
        self.assertIn("rr3", self.model.listOfRateRules)

    def test_model_add__reaction(self):
        from gillespy2 import Reaction
        r4 = Reaction(name="r4", reactants={"s1": 1}, products={"s2": 1}, rate="k1")
        self.model.add(r4)
        self.assertIn("r4", self.model.listOfReactions)

    def test_model_add__species(self):
        from gillespy2 import Species
        s3 = Species(name="s3", initial_value=29)
        self.model.add(s3)
        self.assertIn("s3", self.model.listOfSpecies)

    def test_model_add__timespan(self):
        from gillespy2 import TimeSpan
        tspan = TimeSpan(range(100))
        self.model.add(tspan)
        self.assertEqual(tspan, self.model.tspan)

    def test_model_add__multiple_components__in_order(self):
        import gillespy2

        s1 = gillespy2.Species(name="s1", initial_value=29)
        k1 = gillespy2.Parameter(name="k1", expression=29)
        r1 = gillespy2.Reaction(name="r1", reactants={"s1": 1}, rate="k1")
        rr1 = gillespy2.RateRule(name="rr1", variable="k1", formula="29")
        ar1 = gillespy2.AssignmentRule(name="ar1", variable="s1", formula="29")
        ea = gillespy2.EventAssignment(name="ea", variable="k1", expression="29")
        et = gillespy2.EventTrigger(expression="t > 29")
        e1 = gillespy2.Event(name="e1", trigger=et, assignments=[ea])
        divide = gillespy2.FunctionDefinition(name="divide", function="x / y", args=["x", "y"])
        tspan = gillespy2.TimeSpan(range(100))

        model = gillespy2.Model(name="Test Model")
        model.add([s1, k1, r1, rr1, ar1, e1, divide, tspan])

        self.assertIn("ar1", model.listOfAssignmentRules)
        self.assertIn("e1", model.listOfEvents)
        self.assertIn("divide", model.listOfFunctionDefinitions)
        self.assertIn("k1", model.listOfParameters)
        self.assertIn("rr1", model.listOfRateRules)
        self.assertIn("r1", model.listOfReactions)
        self.assertIn("s1", model.listOfSpecies)
        self.assertEqual(tspan, model.tspan)

    def test_model_add__multiple_components__in_order(self):
        import gillespy2

        s1 = gillespy2.Species(name="s1", initial_value=29)
        k1 = gillespy2.Parameter(name="k1", expression=29)
        r1 = gillespy2.Reaction(name="r1", reactants={"s1": 1}, rate="k1")
        rr1 = gillespy2.RateRule(name="rr1", variable="k1", formula="29")
        ar1 = gillespy2.AssignmentRule(name="ar1", variable="s1", formula="29")
        ea = gillespy2.EventAssignment(name="ea", variable="k1", expression="29")
        et = gillespy2.EventTrigger(expression="t > 29")
        e1 = gillespy2.Event(name="e1", trigger=et, assignments=[ea])
        divide = gillespy2.FunctionDefinition(name="divide", function="x / y", args=["x", "y"])
        tspan = gillespy2.TimeSpan(range(100))

        model = gillespy2.Model(name="Test Model")
        model.add([ar1, divide, e1, k1, s1, r1, rr1, tspan])

        self.assertIn("ar1", model.listOfAssignmentRules)
        self.assertIn("e1", model.listOfEvents)
        self.assertIn("divide", model.listOfFunctionDefinitions)
        self.assertIn("k1", model.listOfParameters)
        self.assertIn("rr1", model.listOfRateRules)
        self.assertIn("r1", model.listOfReactions)
        self.assertIn("s1", model.listOfSpecies)
        self.assertEqual(tspan, model.tspan)

    def test_delete_assignment_rule(self):
        self.model.delete_assignment_rule('rr2')
        self.assertNotIn('rr2', self.model.listOfAssignmentRules)
        self.assertNotIn('rr2', self.model._listOfAssignmentRules)

    def test_delete_event(self):
        self.model.delete_event('e1')
        self.assertNotIn('e1', self.model.listOfEvents)
        self.assertNotIn('e1', self.model._listOfEvents)

    def test_delete_function_definition(self):
        self.model.delete_function_definition('multiply')
        self.assertNotIn('multiply', self.model.listOfFunctionDefinitions)
        self.assertNotIn('multiply', self.model._listOfFunctionDefinitions)

    def test_delete_parameter(self):
        self.model.delete_parameter('k1')
        self.assertNotIn('k1', self.model.listOfParameters)
        self.assertNotIn('k1', self.model._listOfParameters)

    def test_delete_rate_rule(self):
        self.model.delete_rate_rule('rr1')
        self.assertNotIn('rr1', self.model.listOfRateRules)
        self.assertNotIn('rr1', self.model._listOfRateRules)

    def test_delete_reaction(self):
        self.model.delete_reaction('r1')
        self.assertNotIn('r1', self.model.listOfReactions)
        self.assertNotIn('r1', self.model._listOfReactions)

    def test_delete_species(self):
        self.model.delete_species('s1')
        self.assertNotIn('s1', self.model.listOfSpecies)
        self.assertNotIn('s1', self.model._listOfSpecies)
