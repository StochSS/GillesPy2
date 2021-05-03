import unittest
import tempfile
from gillespy2.core.gillespyError import DirectoryError, SimulationError
from example_models import Example
from gillespy2 import SSACSolver, ODECSolver, TauLeapingCSolver


class TestVariableSolvers(unittest.TestCase):
    model = Example()
    solverSSAC = SSACSolver(model, variable=True)
    solverODEC = ODECSolver(model, variable=True)
    solverTAUC = TauLeapingCSolver(model, variable=True)
    solverlist = [solverSSAC, solverODEC, solverTAUC]

    def test_create(self):
        model = Example()
        solverSSAC = SSACSolver(model)
        solverODEC = ODECSolver(model)
        solverTAUC = TauLeapingCSolver(model)

    def test_file_with_directory_name_exists(self):
        with self.assertRaises(DirectoryError):
            temp = tempfile.NamedTemporaryFile()
            model = Example()
            solverSSAC = SSACSolver(model, temp.name)
            solverODEC = ODECSolver(model, temp.name)
            solverTAUC = TauLeapingCSolver(model, temp.name)

    def test_run_example_precompiled(self):
        for solver in self.solverlist:
            results = self.model.run(solver=solver)

    def test_change_species(self):
        initial_value = self.model.listOfSpecies['Sp'].initial_value
        for solver in self.solverlist:
            results = self.model.run(solver=solver, variables={'Sp':3})
        with self.subTest(msg='Test changed species simulation'):
            self.assertEqual(results['Sp'][0], 3)
        with self.subTest(msg='Test changed species model integrity'):
            self.assertEqual(self.model.listOfSpecies['Sp'].initial_value, initial_value)

    def test_change_parameter(self):
        initial_expression = self.model.listOfParameters['k1'].expression
        for solver in self.solverlist:
            print(solver.name)
            results = self.model.run(solver=solver, variables={'k1':0})
            with self.subTest(msg='Test changed parameter simulation'):
                self.assertEqual(results['Sp'][-1], results['Sp'][0])
            with self.subTest(msg='Test changed parameter model integrity'):
                self.assertEqual(self.model.listOfParameters['k1'].expression, initial_expression)

    def test_invalid_variable(self):
        with self.assertRaises(SimulationError):
            for solver in self.solverlist:
                results = self.model.run(solver=solver, variables={'foobar':0})

    def test_run_example(self):
        notPrecompiled = [SSACSolver, ODECSolver, TauLeapingCSolver]
        for solver in notPrecompiled:
            results = self.model.run(solver=solver)


if __name__ == '__main__':
    unittest.main()
