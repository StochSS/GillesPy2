import unittest
import tempfile
from gillespy2.core.gillespyError import DirectoryError, SimulationError
from example_models import Example
from gillespy2 import VariableSSACSolver


class TestVariableSSACSolver(unittest.TestCase):
    def test_create(self):
        model = Example()
        solver = VariableSSACSolver(model)

    def test_file_with_directory_name_exists(self):
        with self.assertRaises(DirectoryError):
            temp = tempfile.NamedTemporaryFile()
            model = Example()
            solver = VariableSSACSolver(model, temp.name)

    def test_run_example_precompiled(self):
        model = Example()
        solver = VariableSSACSolver(model)
        results = model.run(solver=solver)

    def test_change_species(self):
        model = Example()
        initial_value = model.listOfSpecies['Sp'].initial_value
        solver = VariableSSACSolver(model)
        results = model.run(solver=solver, variables={'Sp':3})
        with self.subTest(msg='Test changed species simulation'):
            self.assertEqual(results['Sp'][0], 3)
        with self.subTest(msg='Test changed species model integrity'):
            self.assertEqual(model.listOfSpecies['Sp'].initial_value, initial_value)

    def test_change_parameter(self):
        model = Example()
        initial_expression = model.listOfParameters['k1'].expression
        solver = VariableSSACSolver(model)
        results = model.run(solver=solver, variables={'k1':0})
        with self.subTest(msg='Test changed parameter simulation'):
            self.assertEqual(results['Sp'][-1], results['Sp'][0])
        with self.subTest(msg='Test changed parameter model integrity'):
            self.assertEqual(model.listOfParameters['k1'].expression, initial_expression)

    def test_invalid_variable(self):
        model = Example()
        solver = VariableSSACSolver(model)
        with self.assertRaises(SimulationError):
            results = model.run(solver=solver, variables={'foobar':0})

    def test_run_example(self):
        model = Example()
        results = model.run(solver=VariableSSACSolver)


if __name__ == '__main__':
    unittest.main()
