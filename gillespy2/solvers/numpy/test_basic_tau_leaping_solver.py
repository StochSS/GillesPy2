import unittest
import tempfile
from gillespy2.core.gillespyError import SolverError, DirectoryError, BuildError, ExecutionError
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver


class TestBasicTauLeapingSolver(unittest.TestCase):

    def test_run_example(self):
        model = Example()
        results = model.run(solver=BasicTauLeapingSolver)


if __name__ == '__main__':
    unittest.main()
