from gillespy2.core import Model, Species, Reaction, Parameter, RateRule
import numpy as np


class SimpleHybridModel(Model):
    def __init__(self, parameter_values=None, init_v=1):
        # initialize Model
        Model.__init__(self, name="Simple_Hybrid_Model")

        # Species
        A = Species(name='A', initial_value=0)
        V = Species(name='V', initial_value=init_v)

        self.add_species([A, V])

        # parameters
        rate1 = Parameter(name='rate1', expression=20.0)
        rate2 = Parameter(name='rate2', expression=10.0)
        rate_rule1 = RateRule(V, "cos(t)")
        self.add_parameter([rate1, rate2])
        self.add_rate_rule(rate_rule1)

        # reactions
        r1 = Reaction(name="r1", reactants={}, products={A: 1},
                                propensity_function="rate1 * V")

        r2 = Reaction(name="r2", reactants={A: 1}, products={},
                                rate=rate2)

        self.add_reaction([r1, r2])
        self.timespan(np.linspace(0, 100, 1001))