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

import unittest
import numpy
import gillespy2

class TestNegativeState(unittest.TestCase):
    def setUp(self):

        def create_hybrid_model(include_reaction=False):
            model = gillespy2.Model(name="test_negative_state")
            # Parameters
            k1 = gillespy2.Parameter(name="k1", expression="1e-05")
            model.add_parameter([k1])
            # Variables
            A = gillespy2.Species(name="A", initial_value=10)
            model.add_species([A])
            # Reactions
            if include_reaction:
                model.add_reaction(gillespy2.Reaction(name="r1", reactants={}, 
                    products={'A': 1}, rate="k1"))
            # Rate Rule
            model.add_rate_rule(gillespy2.RateRule(name="rr1", formula="-1", 
                variable="A"))
            # Timespan
            model.timespan(gillespy2.TimeSpan.arange(1,t=20))
            return model

        self.rr_only_model = create_hybrid_model(include_reaction=False)
        self.rr_rxn_model = create_hybrid_model(include_reaction=True)

    def test_hybrid_python_rr_only(self):
        sol = gillespy2.TauHybridSolver(model=self.rr_only_model)
        result = sol.run()
        self.assertLess(result['A'][-1], -9.9)

    def test_hybrid_python_rr_rxn(self):
        sol = gillespy2.TauHybridSolver(model=self.rr_rxn_model)
        with self.assertRaises(gillespy2.SimulationError) as ex:
            result = sol.run()
        self.assertIn('Negative State detected at begining of step. Species involved in reactions can not be negative.', str(ex.exception))


    def test_hybrid_C_rr_only(self):
        sol = gillespy2.TauHybridCSolver(model=self.rr_only_model)
        result = sol.run()
        self.assertLess(result['A'][-1], -9.9)

    def test_hybrid_C_rr_rxn(self):
        sol = gillespy2.TauHybridCSolver(model=self.rr_rxn_model)
        with self.assertRaises(gillespy2.SimulationError) as ex:
            result = sol.run()
        self.assertEqual('Negative State detected at beginning of step. Species involved in reactions can not be negative.', str(ex.exception))


if __name__ == '__main__':
    unittest.main()

