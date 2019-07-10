import unittest
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver


class TestNumPySSASolver(unittest.TestCase):

    def test_run_example(self):
        model = Example()
        results = model.run(solver=NumPySSASolver)

    # Issue #160
    def test_kwargs_throw_error(self):
        solver = NumPySSASolver
        with self.assertRaises(TypeError):
            solver.run(model=Example(), max_steps=0)


if __name__ == '__main__':
    unittest.main()
