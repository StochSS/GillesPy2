import unittest
import tempfile
from gillespy2.core.gillespyError import SolverError, DirectoryError, BuildError, ExecutionError
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver


class TestBasicTauHybridSolver(unittest.TestCase):

    def test_run_example(self):
        model = Example()
        results = model.run(solver=BasicTauHybridSolver)


if __name__ == '__main__':
    unittest.main()
