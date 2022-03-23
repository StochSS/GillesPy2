"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import sys
sys.path.insert(1, "../")
import unittest

from example_models import RobustModel
from gillespy2.core.gillespyError import *

class TestModel(unittest.TestCase):

	def setup(self):
		self.model = RobustModel()

	def test_delete_assignment_rule(self):
        self.model.delete_assignment_rule('k1')
        self.assertNotIn('k1', self.model.listOfAssignmentRules)
        self.assertNotIn('k1', self.model._listOfAssignmentRules)

    def test_delete_event(self):
        self.model.delete_event('k1')
        self.assertNotIn('k1', self.model.listOfEvents)
        self.assertNotIn('k1', self.model._listOfEvents)

    def test_delete_function_definition(self):
        self.model.delete_function_definition('k1')
        self.assertNotIn('k1', self.model.listOfFunctionDefinitions)
        self.assertNotIn('k1', self.model._listOfFunctionDefinitions)

    def test_delete_parameter(self):
        self.model.delete_parameter('k1')
        self.assertNotIn('k1', self.model.listOfParameters)
        self.assertNotIn('k1', self.model._listOfParameters)

    def test_delete_rate_rule(self):
        self.model.delete_rate_rule('k1')
        self.assertNotIn('k1', self.model.listOfRateRules)
        self.assertNotIn('k1', self.model._listOfRateRules)

    def test_delete_reaction(self):
        self.model.delete_reaction('r1')
        self.assertNotIn('r1', self.model.listOfReactions)
        self.assertNotIn('r1', self._model.listOfReactions)

    def test_delete_species(self):
        self.model.delete_species('s1')
        self.assertNotIn('s1', self.model.listOfSpecies)
        self.assertNotIn('s1', self.model._listOfSpecies)
