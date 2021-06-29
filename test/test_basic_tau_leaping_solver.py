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
import numpy as np
from example_models import Example
from gillespy2 import TauLeapingSolver
from gillespy2.core import results

class TestBasicTauLeapingSolver(unittest.TestCase):

    model = Example()
    results = model.run(solver=TauLeapingSolver, show_labels=False, number_of_trajectories=1)
    labels_results = model.run(solver=TauLeapingSolver, show_labels=True, number_of_trajectories=1)

    def test_return_type(self):
        assert(isinstance(self.results, np.ndarray))
        assert(isinstance(self.results[0], np.ndarray))
        assert(isinstance(self.results[0][0], np.ndarray))
        assert(isinstance(self.results[0][0][0], np.float))

    def test_return_type_show_labels(self):
        assert(isinstance(self.labels_results, results))
        assert(isinstance(self.labels_results['Sp'], np.ndarray))
        assert(isinstance(self.labels_results['Sp'][0], np.float))

if __name__ == '__main__':
    unittest.main()
