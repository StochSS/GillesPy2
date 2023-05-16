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
''' Unit tests module for gillespy2.Model. '''
import sys
sys.path.insert(1, "../")
import unittest

from example_models import create_robust_model
import gillespy2

class TestModel(unittest.TestCase):
    ''' Unit tests class for gillespy2.Model. '''
    def setUp(self):
        self.model = create_robust_model()

    def test_model_add__assignment_rule(self):
        ''' Test Model.add with an assignment rule. '''
        ar1 = gillespy2.AssignmentRule(name="ar1", variable="test_species", formula="29")
        self.model.add(ar1)
        self.assertIn("ar1", self.model.listOfAssignmentRules)

    def test_model_add__event(self):
        ''' Test Model.add with an event. '''
        assignment = gillespy2.EventAssignment(variable="k1", expression="29")
        trigger = gillespy2.EventTrigger(expression="t > 29")
        event = gillespy2.Event(name="e2", trigger=trigger, assignments=[assignment])
        self.model.add(event)
        self.assertIn("e2", self.model.listOfEvents)

    def test_model_add__function_definition(self):
        ''' Test Model.add with a function definition. '''
        divide = gillespy2.FunctionDefinition(name="divide", function="x / y", args=["x", "y"])
        self.model.add(divide)
        self.assertIn("divide", self.model.listOfFunctionDefinitions)

    def test_model_add__invalid_component(self):
        ''' Test Model.add with an invalid component. '''
        with self.assertRaises(gillespy2.ModelError):
            self.model.add("29")

    def test_model_add__parameter(self):
        ''' Test Model.add with a parameter. '''
        parameter = gillespy2.Parameter(name="k3", expression=29)
        self.model.add(parameter)
        self.assertIn("k3", self.model.listOfParameters)

    def test_model_add__rate_rule(self):
        ''' Test Model.add with a rate rule. '''
        rr3 = gillespy2.RateRule(name="rr3", variable="test_species", formula="29")
        self.model.add(rr3)
        self.assertIn("rr3", self.model.listOfRateRules)

    def test_model_add__reaction(self):
        ''' Test Model.add with a reaction. '''
        reaction = gillespy2.Reaction(name="r4", reactants={"s1": 1}, products={"s2": 1}, rate="k1")
        self.model.add(reaction)
        self.assertIn("r4", self.model.listOfReactions)

    def test_model_add__species(self):
        ''' Test Model.add with a species. '''
        species = gillespy2.Species(name="s3", initial_value=29)
        self.model.add(species)
        self.assertIn("s3", self.model.listOfSpecies)

    def test_model_add__timespan(self):
        ''' Test Model.add with a timespan. '''
        tspan = gillespy2.TimeSpan(range(100))
        self.model.add(tspan)
        self.assertEqual(tspan, self.model.tspan)

    def test_model_add__multiple_components__in_order(self):
        ''' Test Model.add with multiple components in proper add order. '''
        spec1 = gillespy2.Species(name="s1", initial_value=29, mode="continuous")
        spec2 = gillespy2.Species(name="s2", initial_value=29)
        parameter = gillespy2.Parameter(name="k1", expression=29)
        reaction = gillespy2.Reaction(name="r1", reactants={"s1": 1}, rate="k1")
        rr1 = gillespy2.RateRule(name="rr1", variable="s1", formula="29")
        ar1 = gillespy2.AssignmentRule(name="ar1", variable="s2", formula="29")
        assignment = gillespy2.EventAssignment(variable="k1", expression="29")
        trigger = gillespy2.EventTrigger(expression="t > 29")
        event = gillespy2.Event(name="e1", trigger=trigger, assignments=[assignment])
        divide = gillespy2.FunctionDefinition(name="divide", function="x / y", args=["x", "y"])
        tspan = gillespy2.TimeSpan(range(100))

        model = gillespy2.Model(name="Test Model")
        model.add([spec1, spec2, parameter, reaction, rr1, ar1, event, divide, tspan])

        self.assertIn("ar1", model.listOfAssignmentRules)
        self.assertIn("e1", model.listOfEvents)
        self.assertIn("divide", model.listOfFunctionDefinitions)
        self.assertIn("k1", model.listOfParameters)
        self.assertIn("rr1", model.listOfRateRules)
        self.assertIn("r1", model.listOfReactions)
        self.assertIn("s1", model.listOfSpecies)
        self.assertEqual(tspan, model.tspan)

    def test_model_add__multiple_components__not_in_order(self):
        ''' Test Model.add with multiple components not in proper add order. '''
        spec1 = gillespy2.Species(name="s1", initial_value=29, mode="continuous")
        spec2 = gillespy2.Species(name="s2", initial_value=29)
        parameter = gillespy2.Parameter(name="k1", expression=29)
        reaction = gillespy2.Reaction(name="r1", reactants={"s1": 1}, rate="k1")
        rr1 = gillespy2.RateRule(name="rr1", variable="s1", formula="29")
        ar1 = gillespy2.AssignmentRule(name="ar1", variable="s2", formula="29")
        assignment = gillespy2.EventAssignment(variable="k1", expression="29")
        trigger = gillespy2.EventTrigger(expression="t > 29")
        event = gillespy2.Event(name="e1", trigger=trigger, assignments=[assignment])
        divide = gillespy2.FunctionDefinition(name="divide", function="x / y", args=["x", "y"])
        tspan = gillespy2.TimeSpan(range(100))

        model = gillespy2.Model(name="Test Model")
        model.add([ar1, divide, event, parameter, spec1, reaction, rr1, tspan, spec2])

        self.assertIn("ar1", model.listOfAssignmentRules)
        self.assertIn("e1", model.listOfEvents)
        self.assertIn("divide", model.listOfFunctionDefinitions)
        self.assertIn("k1", model.listOfParameters)
        self.assertIn("rr1", model.listOfRateRules)
        self.assertIn("r1", model.listOfReactions)
        self.assertIn("s1", model.listOfSpecies)
        self.assertEqual(tspan, model.tspan)

    def test_delete_assignment_rule(self):
        ''' Test model.delete_assignment_rule method. '''
        self.model.delete_assignment_rule('rr2')
        self.assertNotIn('rr2', self.model.listOfAssignmentRules)
        self.assertNotIn('rr2', self.model._listOfAssignmentRules) # pylint: disable=protected-access

    def test_delete_event(self):
        ''' Test model.delete_event method. '''
        self.model.delete_event('e1')
        self.assertNotIn('e1', self.model.listOfEvents)
        self.assertNotIn('e1', self.model._listOfEvents) # pylint: disable=protected-access

    def test_delete_function_definition(self):
        ''' Test model.delete_function_definition method. '''
        self.model.delete_function_definition('multiply')
        self.assertNotIn('multiply', self.model.listOfFunctionDefinitions)
        self.assertNotIn('multiply', self.model._listOfFunctionDefinitions) # pylint: disable=protected-access

    def test_delete_parameter(self):
        ''' Test model.delete_parameter method. '''
        self.model.delete_parameter('k1')
        self.assertNotIn('k1', self.model.listOfParameters)
        self.assertNotIn('k1', self.model._listOfParameters) # pylint: disable=protected-access

    def test_delete_rate_rule(self):
        ''' Test model.delete_rate_rule method. '''
        self.model.delete_rate_rule('rr1')
        self.assertNotIn('rr1', self.model.listOfRateRules)
        self.assertNotIn('rr1', self.model._listOfRateRules) # pylint: disable=protected-access

    def test_delete_reaction(self):
        ''' Test model.delete_reaction method. '''
        self.model.delete_reaction('r1')
        self.assertNotIn('r1', self.model.listOfReactions)
        self.assertNotIn('r1', self.model._listOfReactions) # pylint: disable=protected-access

    def test_delete_species(self):
        ''' Test model.delete_species method. '''
        self.model.delete_species('s1')
        self.assertNotIn('s1', self.model.listOfSpecies)
        self.assertNotIn('s1', self.model._listOfSpecies) # pylint: disable=protected-access
