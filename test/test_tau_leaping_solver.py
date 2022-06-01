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

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import unittest
import numpy as np
from example_models import create_decay, create_decay_no_tspan
from gillespy2 import Results, TauLeapingSolver
from gillespy2.core.gillespyError import SimulationError

class TestBasicTauLeapingSolver(unittest.TestCase):

    def setUp(self):
        self.model = create_decay()
        self.model_no_tspan = create_decay_no_tspan()
        self.solver = TauLeapingSolver(model=self.model)
        self.solver_no_tspan = TauLeapingSolver(model=self.model_no_tspan)

    def test_return_type(self):
        results = self.model.run(
            solver=self.solver, number_of_trajectories=1
        ).to_array()
        labels_results = self.model.run(solver=self.solver, number_of_trajectories=1)
        self.assertIsInstance(results, np.ndarray)
        self.assertIsInstance(results[0], np.ndarray)
        self.assertIsInstance(results[0][0], np.ndarray)
        self.assertIsInstance(results[0][0][0], float)

    def test_return_type_show_labels(self):
        results = self.model.run(
            solver=self.solver, number_of_trajectories=1
        ).to_array()
        labels_results = self.model.run(solver=self.solver, number_of_trajectories=1)
        self.assertIsInstance(labels_results, Results)
        self.assertIsInstance(labels_results['Sp'], np.ndarray)
        self.assertIsInstance(labels_results['Sp'][0], float)

    def test_run_example__with_increment_only(self):
        results = self.solver_no_tspan.run(increment=0.2)

    def test_run_example__with_tspan_only(self):
        results = self.solver.run()

    def test_run_example__with_tspan_and_increment(self):
        with self.assertRaises(SimulationError):
            results = self.solver.run(increment=0.2)

if __name__ == '__main__':
    unittest.main()
