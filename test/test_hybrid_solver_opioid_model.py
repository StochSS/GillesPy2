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

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from datetime import time
import unittest
import numpy as np
import io
import gillespy2
from gillespy2.solvers import TauHybridCSolver, TauHybridSolver

class TestHbridSolverOpioidModel(unittest.TestCase):

    def create_opioid(self):
        # Initialize Model
        model = gillespy2.Model(name="Opioid")

        # Define Variables (GillesPy2.Species)
        S = gillespy2.Species(name='Susceptibles', initial_value=1000)
        P = gillespy2.Species(name='Prescribed_Users', initial_value=0)
        A = gillespy2.Species(name='Addicted', initial_value=0)
        R = gillespy2.Species(name='Rehab', initial_value=0)
        Natural_Deaths = gillespy2.Species(name='Natural_Deaths', initial_value=0)
        Addiction_Deaths = gillespy2.Species(name='Addiction_Deaths', initial_value=0)
        
        # Add Variables to Model
        model.add_species([S,P,A,R,Natural_Deaths,Addiction_Deaths])

        # Define Parameters
        alpha = gillespy2.Parameter(name='alpha', expression= 0.15)
        epsilon = gillespy2.Parameter(name='epsilon', expression= 0.8)
        beta_p = gillespy2.Parameter(name='beta_p', expression= 0.00266)
        beta_a = gillespy2.Parameter(name='beta_a', expression= 0.00094)
        gamma = gillespy2.Parameter(name='gamma', expression= 0.00744)
        zeta = gillespy2.Parameter(name='zeta', expression= 0.2)
        delta = gillespy2.Parameter(name='delta', expression= 0.1)
        sigma = gillespy2.Parameter(name='sigma', expression= 0.9)
        mu = gillespy2.Parameter(name='mu', expression= 0.00729)
        mu_prime = gillespy2.Parameter(name='mu_prime', expression= 0.01159)

        # Add Parameters to Model
        model.add_parameter([alpha, epsilon, beta_p, beta_a, gamma, zeta, delta, sigma, mu, mu_prime])

        # Define Reactions
        SP = gillespy2.Reaction(
            name="SP", reactants={'Susceptibles': 1}, products={'Prescribed_Users': 1}, rate='alpha'
        )
        SA_a = gillespy2.Reaction(name="SA_a", reactants={'Susceptibles': 1}, products={'Addicted': 1}, rate='beta_a')
        SA_p = gillespy2.Reaction(name="SA_p", reactants={'Susceptibles': 1}, products={'Addicted': 1}, rate='beta_p')
        mu_S = gillespy2.Reaction(
            name="mu_S", reactants={'Susceptibles': 1}, products={'Susceptibles': 1, 'Natural_Deaths': 1}, rate='mu'
        )
        PA = gillespy2.Reaction(name="PA", reactants={'Prescribed_Users': 1}, products={'Addicted': 1}, rate='gamma')
        PS = gillespy2.Reaction(
            name="PS", reactants={'Prescribed_Users': 1}, products={'Susceptibles': 1}, rate='epsilon'
        )
        AR = gillespy2.Reaction(name="AR", reactants={'Addicted': 1}, products={'Rehab': 1}, rate='zeta')
        RA = gillespy2.Reaction(name="RA", reactants={'Rehab': 1}, products={'Addicted': 1}, rate='delta')
        RS = gillespy2.Reaction(name="RS", reactants={'Rehab': 1}, products={'Susceptibles': 1}, rate='sigma')
        mu_P = gillespy2.Reaction(
            name="mu_P", reactants={'Prescribed_Users': 1},
            products={'Susceptibles': 1, 'Natural_Deaths': 1}, rate='mu'
        )
        mu_R = gillespy2.Reaction(
            name="m_R", reactants={'Rehab': 1}, products={'Susceptibles': 1, 'Natural_Deaths': 1}, rate='mu'
        )
        mu_prime_A = gillespy2.Reaction(
            name="mu_prime_A", reactants={'Addicted': 1},
            products={'Susceptibles': 1, 'Addiction_Deaths': 1}, rate='mu_prime'
        )

        # Add Reactions to Model
        model.add_reaction([SP, PS, SA_a, SA_p, PA, AR, RA, RS, mu_P, mu_R, mu_prime_A, mu_S])
        
        # Define Timespan
        tspan = gillespy2.TimeSpan.linspace(t=10, num_points=11)
        
        # Set Model Timespan
        model.timespan(tspan)
        return model

    def test_run_output(self):
        model = self.create_opioid()
        for solver in [TauHybridSolver, TauHybridCSolver]:
            results = model.run(solver=solver, number_of_trajectories=3, seed=1024)
            for result in results:
                with self.subTest("Processing simulation output for solver {solver.name} for trajectory={n}"):
                    min_s = min(result['Susceptibles'])
                    self.assertGreater(min_s, 500)
                    max_p = max(result['Prescribed_Users'])
                    self.assertLess(max_p, 500)

