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

import unittest
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from profiling import c_profiler

from example_models import Oregonator 
from example_models import VilarOscillator
from example_models import Tyson2StateOscillator 

from gillespy2.solvers.cpp import SSACSolver
from gillespy2.solvers.cpp import ODECSolver
from gillespy2.solvers.cpp import TauLeapingCSolver


class MyTestCase(unittest.TestCase):
    def test_print_profiler_results(self):
        for solver, models in self.solvers.items():
            print(f"=== === === Solver: {solver.name} === === ===")
            for model in models:
                perf_results = c_profiler.run_profiler(model, solver(model=model))
                print(f"  === === Model: {model.name} === ===")
                print(perf_results)
        print(f"=== === === === === === === === === === === === ===")

    def setUp(self) -> None:
        self.solvers = {
            SSACSolver: [
                Tyson2StateOscillator()
            ],
            ODECSolver: [
                Oregonator()
            ],
            TauLeapingCSolver: [
                VilarOscillator()
            ],
        }


if __name__ == '__main__':
    unittest.main()
