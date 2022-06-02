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
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from profiling import c_profiler

from example_models import create_oregonator 
from example_models import create_vilar_oscillator
from example_models import create_tyson_2_state_oscillator 

from gillespy2.solvers.cpp import SSACSolver
from gillespy2.solvers.cpp import ODECSolver
from gillespy2.solvers.cpp import TauLeapingCSolver


class MyTestCase(unittest.TestCase):
    def test_print_profiler_results(self):
        for solver_name, models in self.solvers.items():
            print(f"=== === === \tSolver: {solver_name} \t=== === ===")
            for model in models:
                print(f"[Model: {model.name}]")
                perf_results = c_profiler.run_profiler(model, solver_name)
                print("Call List:")
                for entry, time in list(perf_results.call_list.items())[:10]:
                    print(f"* {time.perf_time:0.1f}ms:\t{entry}")
                print("Execution Time:", f"{perf_results.execution_time:0.1f}ms")
                print("CPU Time:", f"{perf_results.sample_time:0.1f}ms")
        print(f"=== === === === === === === === === === === === ===")

    def setUp(self) -> None:
        self.solvers = {
            SSACSolver.target: [
                create_tyson_2_state_oscillator()
            ],
            ODECSolver.target: [
                create_oregonator()
            ],
        }


if __name__ == '__main__':
    unittest.main()
