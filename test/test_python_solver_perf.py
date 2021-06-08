import pstats
import cProfile
import unittest

from unittest import TestCase

from test.example_models import MichaelisMenten
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
from gillespy2.solvers.numpy.tau_hybrid_solver import TauHybridSolver

class TestPythonSolverPerf(TestCase):
    def setUp(self) -> None:
        self.models = [
            MichaelisMenten()
        ]

        self.solvers = [
            NumPySSASolver(),
            TauHybridSolver()
        ]

        # Generate all possible permutations of self.models and self.solvers.
        self.runnables = [(model, solver) for model in self.models for solver in self.solvers]

    def test_python_solver_perf(self):
        for runnable in self.runnables:
            model = runnable[0]
            solver = runnable[1]

            print(f"Profiling solver: {type(solver)} with model: {type(model)}...")

            profiler = cProfile.Profile()
            profiler.enable()

            solver.run(model=model, number_of_trajectories=100, timeout=100)

            profiler.disable()

            stats = pstats.Stats(profiler).sort_stats("ncalls")
            stats.print_stats()

if __name__ == "__main__":
    unittest.main()
