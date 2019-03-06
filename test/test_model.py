import unittest
from gillespy2.core import Model, Species, Reaction, Parameter
from gillespy2.core.gillespyError import *
import numpy as np


class TestModel(unittest.TestCase):
    def test_duplicate_parameter_names(self):
        model = Model()
        param1 = Parameter('A', expression=0)
        param2 = Parameter('A', expression=0)
        model.add_parameter(param1)
        with self.assertRaises(ModelError):
            model.add_parameter(param2)

    def test_duplicate_species_names(self):
        model = Model()
        species1 = Species('A', initial_value=0)
        species2 = Species('A', initial_value=0)
        model.add_species(species1)
        with self.assertRaises(ModelError):
            model.add_species(species2)

    def test_duplicate_reaction_name(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=0)
        species2 = Species('B', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = Reaction(name="reaction1", reactants={species1: 1}, products={species2: 1}, rate=rate)
        reaction2 = Reaction(name="reaction1", reactants={species2: 1}, products={species1: 1}, rate=rate)
        model.add_reaction(reaction1)
        with self.assertRaises(ModelError):
            model.add_reaction(reaction2)

    def test_species_parameter_name_substrings(self):
        model = Model()
        rate = Parameter(name='rate', expression=1)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=100)
        species2 = Species('AA', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = Reaction(name="reaction1", reactants={species1: 1}, products={species2: 1}, rate=rate)
        model.add_reaction(reaction1)
        number_points = 11
        model.timespan(np.linspace(0, 1, number_points))
        from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
        results = model.run(number_of_trajectories=1, solver=NumPySSASolver, seed=1)[0]
        self.assertTrue(len(results['time']) == number_points)
        self.assertTrue(len(results[species1.name]) == number_points)
        self.assertTrue(len(results[species2.name]) == number_points)
        self.assertGreater(results[species1.name][0], results[species1.name][-1])
        self.assertLess(results[species2.name][0], results[species2.name][-1])
        self.assertEqual(np.sum(results[species1.name]) + np.sum(results[species2.name]), number_points * species1.initial_value)


if __name__ == '__main__':
    unittest.main()
