import unittest
import tempfile
from gillespy2.core.gillespyError import DirectoryError
from example_models import Example
from gillespy2.solvers.cpp.ode_c_solver import ODECSolver


class TestODECSolver(unittest.TestCase):
    def test_create(self):
        model = Example()
        solver = ODECSolver(model)

    def test_file_with_directory_name_exists(self):
        with self.assertRaises(DirectoryError):
            temp = tempfile.NamedTemporaryFile()
            model = Example()
            solver = ODECSolver(model, temp.name)

    def test_run_example_precompiled(self):
        model = Example()
        solver = ODECSolver(model)
        results = model.run(solver=solver)

    # def test_run_example(self):
    #     model = Example()
    #     results = model.run(solver=ODECSolver)


if __name__ == '__main__':
    unittest.main()
