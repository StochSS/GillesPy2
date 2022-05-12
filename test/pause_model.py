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

import sys
sys.path.append('..')
from gillespy2.core import Model, Species, Reaction, Parameter
import numpy as np
from gillespy2 import NumPySSASolver
from gillespy2 import TauLeapingSolver
from gillespy2 import ODESolver
from gillespy2 import TauHybridSolver
import sys
np.set_printoptions(suppress=True)

class Oregonator(Model):

    def __init__(self, parameter_values=None):
        # Superclass initialization
        Model.__init__(self, name="Oregonator")

        # Species
        F = Species(name="F", initial_value=2)
        A = Species(name="A", initial_value=250)
        B = Species(name="B", initial_value=500)
        C = Species(name="C", initial_value=1000)
        P = Species(name="P", initial_value=0)
        self.add_species([F, A, B, C, P])

        # Parameters (rates)
        k1 = Parameter(name="k1", expression=2.0)
        k2 = Parameter(name="k2", expression=0.1)
        k3 = Parameter(name="k3", expression=104)
        k4 = Parameter(name="k4", expression=4e-7)
        k5 = Parameter(name="k5", expression=26.0)
        self.add_parameter([k1, k2, k3, k4, k5])

        # Reactions
        reaction1 = Reaction(name="reaction1",
                             reactants={B: 1, F: 1},
                             products={A: 1, F: 1},
                             rate=k1)
        reaction2 = Reaction(name="reaction2",
                             reactants={A: 1, B: 1},
                             products={P: 1},
                             rate=k2)
        reaction3 = Reaction(name="reaction3",
                             reactants={A: 1, F: 1},
                             products={A: 2, C: 1, F: 1},
                             rate=k3)
        reaction4 = Reaction(name="reaction4",
                             reactants={A: 2},
                             products={P: 1},
                             rate=k4)
        reaction5 = Reaction(name="reaction5",
                             reactants={C: 1, F: 1},
                             products={B: 1, F: 1},
                             rate=k5)
        self.add_reaction([reaction1, reaction2, reaction3, reaction4, reaction5])
        if sys.argv[1] != 'ODESolver':
            self.timespan(np.linspace(0, 5, 501))
        else:
            self.timespan(np.linspace(0, 5, 500001))

def make_solver(which_solver: str, model: Model):
    solver_dict = {solver.name: solver for solver in [
        NumPySSASolver,
        TauLeapingSolver,
        ODESolver,
        TauHybridSolver,
    ]}
    solver_class = solver_dict.get(which_solver)
    return solver_class(model=model)

def run_test(which_solver: str):
    solver = make_solver(which_solver, model=Oregonator())
    results = solver.run()
    print(results.to_array()[0][-1][0])

if __name__ == "__main__":
    run_test(sys.argv[1])
