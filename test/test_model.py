import unittest
import gillespy2
from gillespy2.core import gillespyError
import numpy as np


class TestModel(unittest.TestCase):
    def test_duplicate_parameter_names(self):
        model = gillespy2.Model()
        param1 = gillespy2.Parameter('A', expression=0)
        param2 = gillespy2.Parameter('A', expression=0)
        model.add_parameter(param1)
        with self.assertRaises(gillespyError.ModelError):
            model.add_parameter(param2)

    def test_duplicate_species_names(self):
        model = gillespy2.Model()
        species1 = gillespy2.Species('A', initial_value=0)
        species2 = gillespy2.Species('A', initial_value=0)
        model.add_species(species1)
        with self.assertRaises(gillespyError.ModelError):
            model.add_species(species2)

    def test_duplicate_reaction_name(self):
        model = gillespy2.Model()
        rate = gillespy2.Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = gillespy2.Species('A', initial_value=0)
        species2 = gillespy2.Species('B', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = gillespy2.Reaction(name="reaction1", reactants={species1: 1}, products={species2: 1}, rate=rate)
        reaction2 = gillespy2.Reaction(name="reaction1", reactants={species2: 1}, products={species1: 1}, rate=rate)
        model.add_reaction(reaction1)
        with self.assertRaises(gillespyError.ModelError):
            model.add_reaction(reaction2)

    def test_species_parameter_name_substrings(self):
        model = gillespy2.Model()
        rate = gillespy2.Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = gillespy2.Species('A', initial_value=1)
        species2 = gillespy2.Species('AA', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = gillespy2.Reaction(name="reaction1", reactants={species1: 1}, products={species2: 1}, rate=rate)
        reaction2 = gillespy2.Reaction(name="reaction2", reactants={species2: 1}, products={species1: 1}, rate=rate)
        model.add_reaction([reaction1, reaction2])
        number_points = 11
        model.timespan(np.linspace(0, 10, number_points))
        results = model.run(number_of_trajectories=1)[0]
        self.assertTrue(len(results['time']) == number_points)
        self.assertTrue(len(results[species1.name]) == number_points)
        self.assertTrue(len(results[species2.name]) == number_points)
        self.assertGreater(np.sum(results[species1.name]), 0)
        self.assertGreater(np.sum(results[species2.name]), 0)
        self.assertEqual(np.sum(results[species1.name]) + np.sum(results[species2.name]), number_points)





