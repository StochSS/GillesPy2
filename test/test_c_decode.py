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
from gillespy2.solvers import ODECSolver, SSACSolver, TauLeapingCSolver, TauHybridCSolver
from example_models import create_degradation
from gillespy2.solvers.cpp.c_decoder import BasicSimDecoder, IterativeSimDecoder

class TestCDecode(unittest.TestCase):
    decoders = [BasicSimDecoder, IterativeSimDecoder]

    def setUp(self):
        self.model = create_degradation()
        tspan = gillespy2.TimeSpan.linspace(t=300, num_points=301)
        self.model.timespan(tspan)
        self.solvers = [
            # Solver and List of trajectory numbers to try
            (ODECSolver(model=self.model), [1]),
            (SSACSolver(model=self.model), [1, 3]),
            (TauLeapingCSolver(model=self.model), [1, 3]),
            (TauHybridCSolver(model=self.model), [1, 3]),
        ]


    def test_c_decoder(self):
        """
        Ensures that, given a certain input from stdout, the results will be
            properly formatted by the sim decoders.
        """
        example_input = """
        0,1.0,2.0,3.0,1,2.0,4.0,6.0,2,3.0,6.0,9.0,3,4.0,8.0,12.0,
        0,1.0,2.0,3.0,1,2.0,4.0,6.0,2,3.0,6.0,9.0,3,4.0,8.0,12.0,3
        """.strip()
        expected_result = np.array([
            [ [0, 1.0,2.0,3.0], [1, 2.0,4.0,6.0], [2, 3.0,6.0,9.0], [3, 4.0,8.0,12.0] ],
            [ [0, 1.0,2.0,3.0], [1, 2.0,4.0,6.0], [2, 3.0,6.0,9.0], [3, 4.0,8.0,12.0] ]
        ])

        for decoder in TestCDecode.decoders:
            mock_stdout = io.BytesIO(bytes(example_input, "ascii"))
            mock_stdout = io.BufferedReader(mock_stdout)
            with self.subTest("Processing mock data output with C decoders", decoder=decoder):
                reader = decoder(np.zeros((2, 4, 4)))
                reader.read(mock_stdout)
                result, time_stopped = reader.get_output()

                self.assertEqual(time_stopped, 3, f"Stop time returned by decoder is incorrect (expected {3}, got {time_stopped})")
                # Test both the passed-in trajectory list and the returned trajectories.
                self.assertTrue(np.all(result == expected_result), "Trajectory data returned by decoder is incorrect")
