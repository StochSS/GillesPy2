import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import unittest
from unittest import TestCase

from profiling import python_profiler

from example_models import Oregonator
from example_models import VilarOscillator
from example_models import Tyson2StateOscillator

from gillespy2.solvers.numpy.ode_solver import ODESolver
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
from gillespy2.solvers.numpy.tau_leaping_solver import TauLeapingSolver


class TestPythonSolverPerf(TestCase):
    def setUp(self) -> None:
        self.solvers = {
            NumPySSASolver: [
                Tyson2StateOscillator()
            ],
            ODESolver: [
                Oregonator()
            ],
            TauLeapingSolver: [
                VilarOscillator()
            ],
        }

    def test_python_solver_perf(self):
        for solver, models in self.solvers.items():
            print(f"=== === === {solver.name} === === ===")

            for model in models:
                print(f"{model.name}:")

                perf_data = python_profiler.run_profiler(model, solver)
                print(perf_data)

            print()

        print()

if __name__ == "__main__":
    unittest.main()
