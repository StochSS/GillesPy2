import unittest

from unittest import TestCase

from test.example_models import MichaelisMenten
from test.perf import python_profiler
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

            python_profiler.run_profiler(model, solver)

if __name__ == "__main__":
    unittest.main()
