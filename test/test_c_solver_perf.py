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
