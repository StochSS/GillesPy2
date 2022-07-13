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
from example_models import *
from gillespy2.core.gillespyError import *
from gillespy2 import ODESolver
#from TestBattery import *

class TestExampleModels(unittest.TestCase):

    def test_trichloroethylene_example(self):
        trichloroethylene_model = create_trichloroethylene()
        results = trichloroethylene_model.run()

    def test_lacOperon_example(self):
        lacOperon_model = create_lac_operon()
        solver = ODESolver(model=lacOperon_model)
        results = lacOperon_model.run(solver=solver)

    def test_schlogl_example(self):
        schlogl_model = create_schlogl()
        results = schlogl_model.run()

    def test_michaelisMenten_example(self):
        michaelisMenten_model = create_michaelis_menten()
        results = michaelisMenten_model.run()

    def test_toggleSwitch_example(self):
        toggleSwitch_model = create_toggle_switch()
        results = toggleSwitch_model.run()

    def test_example_example(self):
        example_model = create_decay()
        results = example_model.run()

    def test_tyson2StateOscillator_example(self):
        tyson2StateOscillator_model = create_tyson_2_state_oscillator()
        results = tyson2StateOscillator_model.run()
    
    #def test_test_battery(self):
        #timing_battery(100, 0.1)

if __name__ == '__main__':
    unittest.main()
