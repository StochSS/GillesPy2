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

from datetime import time
import unittest
import numpy as np
import io
import gillespy2
from gillespy2.solvers import ODESolver, ODECSolver, NumPySSASolver, SSACSolver, TauLeapingCSolver, TauLeapingSolver, TauHybridSolver, TauHybridCSolver
from example_models import create_degradation

class TestRunOutput(unittest.TestCase):

    def setUp(self):
        self.model = create_degradation()
        tspan = gillespy2.TimeSpan.linspace(t=300, num_points=301)
        self.model.timespan(tspan)
        self.solvers = [
            # Solver and List of trajectory numbers to try
            (ODECSolver(model=self.model), [1]),
            (SSACSolver(model=self.model), [1, 3]),
            (TauLeapingCSolver(model=self.model), [1, 3]),
            (TauHybridCSolver(model=self.model), [1, 3]),
            (ODESolver(model=self.model), [1]),
            (NumPySSASolver(model=self.model), [1, 3]),
            (TauLeapingSolver(model=self.model), [1, 3]),
            (TauHybridSolver(model=self.model), [1, 3]),
        ]

    def test_run_output(self):
        expected_tspan = self.model.tspan
        expected_init = self.model.get_species("A").initial_value
        for solver, trajectory_counts in self.solvers:
            for number_of_trajectories in trajectory_counts:
                with self.subTest("Processing simulation output for each solver with different trajectory sizes",
                                  number_of_trajectories=number_of_trajectories,
                                  solver=solver):
                    results = self.model.run(solver=solver, number_of_trajectories=number_of_trajectories, seed=1024)
                    for n,result in enumerate(results):
                        tspan, first_value, last_value = result["time"], result["A"][0], result["A"][-1]
                        result_diff = np.concatenate([np.array([first_value]), result["A"][1:]])
                        result_diff = result["A"] - result_diff
                        self.assertTrue(np.allclose(tspan, expected_tspan), msg="Simulation output contains unexpected timeline values"
                                                                        f"\n  Received timeline: {tspan}\n  Expected timeline: {expected_tspan}, trajectory_number={n}")
                        self.assertEqual(first_value, expected_init, msg=f"Simulation output begins with unexpected value: {first_value}, trajectory_number={n}")
                        self.assertAlmostEqual(last_value, 0, places=3, msg=f"Simulation output converges on an unexpectedly large value: {last_value}, trajectory_number={n}")

