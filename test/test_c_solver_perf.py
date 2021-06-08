import unittest
from test.example_models import Tyson2StateOscillator, VilarOscillator, Oregonator
from gillespy2.solvers.cpp import SSACSolver, ODECSolver, TauLeapingCSolver
from test.perf.gprof import run_profiler


class MyTestCase(unittest.TestCase):
    def test_print_profiler_results(self):
        for solver, models in self.solvers.items():
            print(f"=== === === Solver: {solver.name} === === ===")
            for model in models:
                perf_results = run_profiler(model, solver(model=model))
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
