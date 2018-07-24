import unittest
import tempfile
from gillespy2 import SSACMultiSolver, InvalidAlphaError, InvalidProcessesError
from gillespy2.example_models import Example


class TestSSACMultiSolver(unittest.TestCase):
    def test_create(self):
        model = Example()
        solver = SSACMultiSolver(model=model, alpha=1, win_py_native=True)

    def test_alpha_greater_than_1(self):
        with self.assertRaises(InvalidAlphaError):
            temp = tempfile.NamedTemporaryFile()
            model = Example()
            solver = SSACMultiSolver(model=model, alpha=1.1, win_py_native=True)

    def test_alpha_less_than_equal_to_zero(self):
        with self.assertRaises(InvalidAlphaError):
            temp = tempfile.NamedTemporaryFile()
            model = Example()
            solver = SSACMultiSolver(model=model, alpha=0, win_py_native=True)

    def test_processes_less_than_equal_to_zero(self):
            temp = tempfile.NamedTemporaryFile()
            model = Example()
            solver = SSACMultiSolver(model=model, number_of_processes=0, win_py_native=True)

    def test_processes_greater_than_100(self):
        with self.assertRaises(InvalidProcessesError):
            temp = tempfile.NamedTemporaryFile()
            model = Example()
            solver = SSACMultiSolver(model=model, number_of_processes=101, win_py_native=True)

    def test_run_example_precompiled(self):
        model = Example()
        solver = SSACMultiSolver(model=model, alpha=1, win_py_native=True)
        results = model.run(solver=solver)

    def test_run_example(self):
        model = Example()
        results = model.run(solver=SSACMultiSolver(model=model, alpha=1, win_py_native=True))

        """Tests to make:
        invalid input alpha :: Create AlphaError
        1-100 processes in loop :: Create InvalidNumberProcessesError
        0 and - processes
        processes > 100
        output_directory && delete_directory specified
        """


if __name__ == '__main__':
    unittest.main()