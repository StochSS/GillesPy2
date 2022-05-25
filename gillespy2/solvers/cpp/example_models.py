# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from gillespy2.core import Model, Species, Reaction, Parameter
import numpy as np


def build_model(parameter_values=None):
    """
    This is a simple example for mass-action degradation of species S.
    """
    # Initialize the model.
    model = Model(name="Example")
    # Species
    S = Species(name='Sp', initial_value=100)
    model.add_species([S])
    # Parameters
    k1 = Parameter(name='k1', expression=3.0)
    model.add_parameter([k1])
    # Reactions
    rxn1 = Reaction(name='S degradation', reactants={S: 1}, products={}, rate=k1)
    model.add_reaction([rxn1])
    model.timespan(np.linspace(0, 20, 101))
    return model


__all__ = ['Trichloroethylene', 'LacOperon', 'Schlogl', 'MichaelisMenten',
           'ToggleSwitch', 'build_model', 'Tyson2StateOscillator', 'Oregonator']
