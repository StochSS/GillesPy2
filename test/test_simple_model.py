import unittest
import sys, os
sys.path.append(os.path.abspath(os.getcwd()))
import numpy, math, gillespy2
from gillespy2.basic_tau_hybrid_solver import BasicTauHybridSolver
import matplotlib.pyplot as plt


class SimpleHybridModel(gillespy2.Model):
     def __init__(self, parameter_values=None):
            #initialize Model
            gillespy2.Model.__init__(self, name="Simple_Hybrid_Model")

            # Species
            A = gillespy2.Species(name='A', initial_value=0)
            B = gillespy2.Species(name='B', initial_value=0)
            self.add_species([A,B])

            # Parameters
            k1 = gillespy2.Parameter(name='k1', expression=1)
            k2 = gillespy2.Parameter(name='k2', expression=10)
            self.add_parameter([k1,k2])

            # Rate Rule
            rate_rule = gillespy2.RateRule(B, "cos(t)")
            self.add_rate_rule(rate_rule)

            # Reactions
            r1 = gillespy2.Reaction(name='r1', reactants={A:1}, products={}, propensity_function="k1*B")
            r2 = gillespy2.Reaction(name='r2', reactants={A:1}, products={}, rate=k2)
            self.add_reaction([r1,r2])

            self.timespan(numpy.linspace(0,1,11))

class TestSimpleModel(unittest.TestCase):
    def setUp(self):
        self.model = SimpleHybridModel()
    def test_model_creation(self):
        name = self.model.name
        self.assertEqual(name, "Simple_Hybrid_Model", msg="Unexpected value: {}".format(name))

    def test_model_has_species(self):
        species = self.model.get_species('A')
        self.assertIsInstance(species, gillespy2.Species, msg='{0} has incorrect type'.format(species))

    def test_model_has_parameters(self):
        parameters = self.model.get_all_parameters()
        self.assertIsInstance(parameters['k1'], gillespy2.Parameter, msg='{0} has incorrect type'.format(parameters))
        self.assertIsInstance(parameters['k2'], gillespy2.Parameter, msg='{0} has incorrect type'.format(parameters))

    def test_model_parameters_correct(self):
        parameters = self.model.get_all_parameters()
        self.assertEqual(parameters['k1'].expression, '1', msg='Has incorrect expression')
        self.assertEqual(parameters['k2'].expression, '10', msg='Has incorrect expression')

    def test_model_has_rate_rules(self):
        rate_rules = self.model.listOfRateRules
        self.assertEqual(rate_rules['B'].species.name, 'B', msg='Has incorrect species')
        self.assertEqual(rate_rules['B'].expression, 'cos(t)', msg='{0} has incorrect type'.format(rate_rules))

    def test_model_has_reactions(self):
        reactions = self.model.get_all_reactions()
        self.assertIsInstance(reactions['r1'], gillespy2.Reaction, msg='{0} has incorrect type'.format(reactions))
        self.assertIsInstance(reactions['r2'], gillespy2.Reaction, msg='{0} has incorrect type'.format(reactions))

    def test_model_has_reactions_correct(self):
        reactions = self.model.get_all_reactions()

        species_A = self.model.get_species('A')
        reactants_r1 = reactions['r1'].reactants
        species_key = list(reactants_r1)[0]

        self.assertEqual(reactions['r1'].name, 'r1', msg='Has incorrect expression')
        self.assertEqual(reactions['r1'].propensity_function, 'k1*B', msg='Has incorrect expression')

        self.assertEqual(reactions['r2'].name, 'r2', msg='Has incorrect expression')
        self.assertEqual(reactions['r2'].marate.expression, '10', msg='Has incorrect expression')

        self.assertIsInstance(species_key, gillespy2.Species, msg='Has incorrect type')
        self.assertEqual(reactants_r1[species_A], 1, msg='Has incorrect number of reactants')
    
    def test_model_has_timespan_correct(self):
        timespan = self.model.tspan
        self.assertCountEqual(timespan, numpy.linspace(0,1,11), msg='Has incorrect timespan')

if __name__ == '__main__':
    unittest.main()