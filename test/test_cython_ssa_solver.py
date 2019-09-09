import unittest
from gillespy2.example_models import Example

from gillespy2.solvers.cython import can_use_cython
if can_use_cython:
    from gillespy2.solvers.cython.cython_ssa_solver import CythonSSASolver


class TestCythonSSASolver(unittest.TestCase):

    def test_run_example(self):
        if can_use_cython:
            model = Example()
            results = model.run(solver=CythonSSASolver)


if __name__ == '__main__':
    unittest.main()
