import unittest
import numpy as np
import gillespy2
from gillespy2.core.gillespyError import *
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver


class TestBasicTauHybridSolver(unittest.TestCase):
    model = Example()
    
    def test_add_rate_rule(self):
        species = gillespy2.Species('test_species', initial_value=1, mode='continuous')
        rule = gillespy2.RateRule(species, 'cos(t)')
        self.model.add_species([species])
        self.model.add_rate_rule([rule])
        self.model.run(solver=BasicTauHybridSolver)

    def test_add_rate_rule_dict(self):
        species2 = gillespy2.Species('test_species2', initial_value=2, mode='continuous')
        species3 = gillespy2.Species('test_species3', initial_value=3, mode='continuous')
        rule2 = gillespy2.RateRule(species2, 'cos(t)')
        rule3 = gillespy2.RateRule(variable=species3, formula='sin(t)')
        rate_rule_dict = {'rule2': rule2, 'rule3': rule3}
        self.model.add_species([species2, species3])
        with self.assertRaises(ParameterError):
            self.model.add_rate_rule(rate_rule_dict)
        
    def test_add_bad_species_rate_rule_dict(self):
        species2 = gillespy2.Species('test_species2', initial_value=2, mode='continuous')
        rule = gillespy2.RateRule(formula='sin(t)')
        with self.assertRaises(ModelError):
            self.model.add_rate_rule(rule)

    def test_add_bad_expression_rate_rule_dict(self):
        species2 = gillespy2.Species('test_species2', initial_value=2, mode='continuous')
        rule = gillespy2.RateRule(variable=species2, formula='')
        with self.assertRaises(ModelError):
            self.model.add_rate_rule(rule)


if __name__ == '__main__':
    unittest.main()
