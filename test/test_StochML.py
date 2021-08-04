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
from example_models import Example
from gillespy2.core import Model, StochMLDocument
from gillespy2.core.gillespyError import *
from gillespy2 import ODESolver
from gillespy2 import NumPySSASolver


class TestStochML(unittest.TestCase):

    def test_StochML_from_and_to_model(self):
        model = Example()
        stochml = StochMLDocument.from_model(model)
        stochml_model = stochml.to_model('model')
        stochml_model.run(solver=ODESolver)
        stochml_model.run(solver=NumPySSASolver)



if __name__ == '__main__':
    unittest.main()
