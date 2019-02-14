import unittest
from gillespy2.example_models import Example
from gillespy2.solvers.python.basic_ssa_solver import BasicSSASolver


class TestBasicSSASolver(unittest.TestCase):

    def test_run_example(self):
        model = Example()
        results = model.run(solver=BasicSSASolver)


if __name__ == '__main__':
    unittest.main()
