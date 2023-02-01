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
import gillespy2
from example_models import create_telegraph_model
from gillespy2.core import Model, Species, Reaction, Parameter
from gillespy2.core.gillespyError import *
import os
import numpy as np


class TestTelegraphModel(unittest.TestCase):

    def test_telegraph_model(self):
        model = create_telegraph_model()
        solvers = [
            gillespy2.solvers.CLESolver,
            gillespy2.solvers.NumPySSASolver,
            gillespy2.solvers.ODECSolver,
            gillespy2.solvers.ODESolver,
            gillespy2.solvers.SSACSolver,
            gillespy2.solvers.TauHybridCSolver,
            gillespy2.solvers.TauHybridSolver,
            gillespy2.solvers.TauLeapingCSolver,
            gillespy2.solvers.TauLeapingSolver,
        ]
        for solver in solvers:
            with self.subTest(solver=solver.name):
                results = model.run(solver=solver)
        
        a = np.sum(results['ON'])/len(results['ON'])
        b = np.sum(results['OFF'])/len(results['OFF'])
        self.assertTrue( a < b) 



if __name__ == '__main__':
    unittest.main()
