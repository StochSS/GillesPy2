import unittest
import tempfile
from gillespy2.core.gillespyError import DirectoryError
from gillespy2.example_models import Example
from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver


class TestSSACSolver(unittest.TestCase):
    def test_create(self):
        model = Example()
        solver = SSACSolver(model)

    def test_file_with_directory_name_exists(self):
        with self.assertRaises(DirectoryError):
            temp = tempfile.NamedTemporaryFile()
            model = Example()
            solver = SSACSolver(model, temp.name)

    def test_run_example_precompiled(self):
        model = Example()
        solver = SSACSolver(model)
        results = model.run(solver=solver)

    def test_run_example(self):
        model = Example()
        results = model.run(solver=SSACSolver)


if __name__ == '__main__':
    unittest.main()
