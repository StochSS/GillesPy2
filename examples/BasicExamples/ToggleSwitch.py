from gillespy2.core import Model, Species, Reaction, Parameter
import numpy as np


class ToggleSwitch(Model):
    """
    Gardner et al. Nature (1999)Construction of a genetic toggle switch in Escherichia coli
    (Transcription from
    """

    def __init__(self, parameter_values=None):
        # Initialize the model.
        Model.__init__(self, name="Toggle_Switch")
        # Species
        A = Species(name='A', initial_value=10.0)
        B = Species(name='B', initial_value=10.0)
        self.add_species([A, B])
        # Parameters
        alpha1 = Parameter(name='alpha1', expression=10.0)
        alpha2 = Parameter(name='alpha2', expression=10.0)
        beta = Parameter(name='beta', expression=2.0)
        gamma = Parameter(name='gamma', expression=2.0)
        mu = Parameter(name='mu', expression=1.0)
        self.add_parameter([alpha1, alpha2, beta, gamma, mu])
        # Reactions
        cu = Reaction(name="r1", reactants={}, products={A: 1},
                      rate=alpha1.value * (1 + B.initial_value ** float(beta.expression)))
        cv = Reaction(name="r2", reactants={}, products={B: 1},
                      rate=alpha2.value(1 + A.initial_value ** float(gamma.expression)))
        du = Reaction(name="r3", reactants={A: 1}, products={}, rate=mu)
        dv = Reaction(name="r4", reactants={B: 1}, products={}, rate=mu)
        self.add_reaction([cu, cv, du, dv])
        self.timespan(np.linspace(0, 250, 251))

