import unittest
import tempfile
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))
from gillespy2.core.gillespyError import SolverError, DirectoryError, BuildError, ExecutionError
from gillespy2.example_models import Example
from gillespy2.solvers.cython.cython_ssa_solver import CythonSSASolver


class TestCythonSSASolver(unittest.TestCase):

    def test_run_example(self):
        model = Example()
        results = model.run(solver=CythonSSASolver)


if __name__ == '__main__':
    unittest.main()
