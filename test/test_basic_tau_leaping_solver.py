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
import numpy as np
from example_models import Example, ExampleNoTspan
from gillespy2 import TauLeapingSolver
from gillespy2.core.gillespyError import SimulationError
from gillespy2.core import results

class TestBasicTauLeapingSolver(unittest.TestCase):

    def setUp(self):
        self.model = Example()
        self.model_no_tspan = ExampleNoTspan()
        self.solver = TauLeapingSolver(model=self.model)
        self.solver_no_tspan = TauLeapingSolver(model=self.model_no_tspan)

    def test_return_type(self):
        results = self.model.run(solver=self.solver, show_labels=False, number_of_trajectories=1)
        labels_results = self.model.run(solver=self.solver, show_labels=True, number_of_trajectories=1)
        assert(isinstance(results, np.ndarray))
        assert(isinstance(results[0], np.ndarray))
        assert(isinstance(results[0][0], np.ndarray))
        assert(isinstance(results[0][0][0], np.float))

    def test_return_type_show_labels(self):
        results = self.model.run(solver=self.solver, show_labels=False, number_of_trajectories=1)
        labels_results = self.model.run(solver=self.solver, show_labels=True, number_of_trajectories=1)
        assert(isinstance(labels_results, results))
        assert(isinstance(labels_results['Sp'], np.ndarray))
        assert(isinstance(labels_results['Sp'][0], np.float))

    def test_run_example__with_increment_only(self):
        results = self.solver_no_tspan.run(increment=0.2)

    def test_run_example__with_tspan_only(self):
        model = Example()
        results = TauLeapingSolver.run(model=model)

    def test_run_example__with_tspan_and_increment(self):
        with self.assertRaises(SimulationError):
            results = self.solver.run(increment=0.2)

if __name__ == '__main__':
    unittest.main()
