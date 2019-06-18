import unittest
import numpy as np
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver


class TestBasicODESolver(unittest.TestCase):
    number_of_trajectories = 10

    def test_run_example(self):
        model = Example()
        for i in [1, self.number_of_trajectories]:
            for label in [True, False]:
                with self.subTest(number_of_trajectories=i, show_labels=label):
                    if i > 1:
                        with self.assertWarns(Warning):
                            results = model.run(solver=BasicODESolver, show_labels=label, number_of_trajectories=i)
                    else:
                        results = model.run(solver=BasicODESolver, show_labels=label, number_of_trajectories=i)
                    self.assertEqual(len(results), i)
                    self.assertTrue(all([isinstance(result, dict if label else np.ndarray) for result in results]))
                    if label:
                        result = results[0]
                        for species in model.listOfSpecies.keys():
                            self.assertIn(species, result.keys())
                            self.assertIsInstance(result[species], np.ndarray)
                            self.assertListEqual(list(result[species].shape), list(model.tspan.shape))
                    else:
                        self.assertIsInstance(results, np.ndarray)
                        self.assertListEqual(list(results.shape), [i, len(model.tspan), len(model.listOfSpecies.keys())+1])
                    for result in results[1:]:
                        if label:
                            self.assertDictEqual(results[0], result)
                        else:
                            self.assertTrue(np.array_equal(results[0], result))


if __name__ == '__main__':
    unittest.main()
