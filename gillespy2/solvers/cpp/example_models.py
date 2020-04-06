from gillespy2.core import Model, Species, Reaction, Parameter
import numpy as np


class Example(Model):
    """
    This is a simple example for mass-action degradation of species S.
    """

    def __init__(self, parameter_values=None):
        # Initialize the model.
        Model.__init__(self, name="Example")
        # Species
        S = Species(name='Sp', initial_value=100)
        self.add_species([S])
        # Parameters
        k1 = Parameter(name='k1', expression=3.0)
        self.add_parameter([k1])
        # Reactions
        rxn1 = Reaction(name='S degradation', reactants={S: 1}, products={}, rate=k1)
        self.add_reaction([rxn1])
        self.timespan(np.linspace(0, 20, 101))


__all__ = ['Trichloroethylene', 'LacOperon', 'Schlogl', 'MichaelisMenten',
           'ToggleSwitch', 'Example', 'Tyson2StateOscillator', 'Oregonator']
