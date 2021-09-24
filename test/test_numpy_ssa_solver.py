"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import unittest
from gillespy2.core.gillespyError import SimulationError
from example_models import Example, ExampleNoTspan
from gillespy2 import NumPySSASolver


class TestNumPySSASolver(unittest.TestCase):
    model = Example()

    def test_run_example__with_increment_only(self):
        model = ExampleNoTspan()
        results = NumPySSASolver.run(model=model, increment=0.2)

    def test_run_example__with_tspan_only(self):
        model = Example()
        results = NumPySSASolver.run(model=model)

    def test_run_example__with_tspan_and_increment(self):
        with self.assertRaises(SimulationError):
            model = Example()
            results = NumPySSASolver.run(model=model, increment=0.2)
    

if __name__ == '__main__':
    unittest.main()
