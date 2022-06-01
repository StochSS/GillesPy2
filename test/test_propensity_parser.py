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

import unittest
from gillespy2.solvers.utilities.solverutils import dependency_grapher
from gillespy2.core import Reaction, Species


class TestPropensityFunctions(unittest.TestCase):
    def test_dependency_graphing(self):
        from example_models import create_toggle_switch, create_michaelis_menten
        model = create_toggle_switch()
        dependencies = dependency_grapher(model, list(model.listOfReactions.keys()))
        correct_graph = {'cu': {'dependencies': ['cv', 'dv']}, 'cv': {'dependencies': ['cu', 'du']},
                         'du': {'dependencies': ['cu']}, 'dv': {'dependencies': ['cv']}}
        self.assertEqual(correct_graph, dependencies)

        model = create_michaelis_menten()
        dependencies = dependency_grapher(model, list(model.listOfReactions.keys()))
        correct_graph = {'r1': {'dependencies': ['r2', 'r3']}, 'r2': {'dependencies': ['r1', 'r3']},
                         'r3': {'dependencies': ['r1', 'r2']}}
        self.assertEqual(correct_graph, dependencies)
