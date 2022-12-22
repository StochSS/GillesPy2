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
from gillespy2.solvers import TauHybridCSolver, TauHybridSolver
from example_models import create_opioid

class TestHbridSolverOpioidModel(unittest.TestCase):

    def test_run_output(self):
        model = create_opioid()
        for solver in [TauHybridSolver, TauHybridCSolver]:
            results = model.run(solver=solver, number_of_trajectories=3, seed=1024)
            for result in results:
                with self.subTest("Processing simulation output for solver {solver.name} for trajectory={n}"):
                    min_s = min(result['Susceptibles'])
                    self.assertGreater(min_s, 500)
                    max_p = max(result['Prescribed_Users'])
                    self.assertLess(max_p, 500)

