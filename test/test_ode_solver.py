import unittest
import numpy as np
import gillespy2
from example_models import Example
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver


class TestBasicODESolver(unittest.TestCase):
    number_of_trajectories = 10

    def test_run_example(self):
        model = Example()
        for i in [1, self.number_of_trajectories]:
            for label in [True, False]:
                with self.subTest(number_of_trajectories=i, show_labels=label):
                    if i > 1:
                        with self.assertLogs(level='WARN'):
                            results = model.run(solver=BasicODESolver, show_labels=label, number_of_trajectories=i)
                        self.assertEqual(len(results), i)
                    else:
                        results = model.run(solver=BasicODESolver, show_labels=label, number_of_trajectories=i)

                    if i > 1 or not label:
                        for result in results[1:]:
                            if label:
                                self.assertEqual(results[0], result)

                            else:
                                self.assertTrue(np.array_equal(results[0], result))

    def test_stoich2(self):
        class StoichTestModel(gillespy2.Model):
            def __init__(self, parameter_values=None):
                gillespy2.Model.__init__(self, name='StochTest1')
                A = gillespy2.Species(name='A', initial_value=10)
                B = gillespy2.Species(name='B', initial_value=0)
                self.add_species([A, B])
                k = gillespy2.Parameter(name='k', expression=10)
                self.add_parameter([k])
                r = gillespy2.Reaction(name='r', reactants={A: 2}, products={B:1}, rate=k)
                self.add_reaction([r])
                self.timespan(np.linspace(0, 100, 101))
        model = StoichTestModel()
        result = model.run(solver=BasicODESolver)
        self.assertAlmostEqual(result['B'][-1], 5, places=3)



if __name__ == '__main__':
    unittest.main()
