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
from unittest import TestCase

from profiling import python_profiler

from example_models import create_oregonator
from example_models import create_vilar_oscillator
from example_models import create_tyson_2_state_oscillator

from gillespy2.solvers.numpy.ode_solver import ODESolver
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
from gillespy2.solvers.numpy.tau_leaping_solver import TauLeapingSolver


class TestPythonSolverPerf(TestCase):
    def setUp(self) -> None:
        self.solvers = {
            NumPySSASolver: [
                create_tyson_2_state_oscillator()
            ],
            ODESolver: [
                create_oregonator()
            ],
            TauLeapingSolver: [
                create_vilar_oscillator()
            ],
        }

    def test_python_solver_perf(self):
        for solver, models in self.solvers.items():
            for model in models:
                with self.subTest(model=model.name, solver=solver.name):
                    perf_data = python_profiler.run_profiler(model, solver)
                    print(perf_data)

if __name__ == "__main__":
    unittest.main()
