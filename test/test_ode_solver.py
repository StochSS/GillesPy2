import unittest
import tempfile
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))
from gillespy2.core.gillespyError import SolverError, DirectoryError, BuildError, ExecutionError
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver


class TestBasicODESolver(unittest.TestCase):

    def test_run_example(self):
        model = Example()
        results = model.run(solver=BasicODESolver)


if __name__ == '__main__':
    unittest.main()
