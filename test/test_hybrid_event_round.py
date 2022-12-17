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
from gillespy2 import Model, Species, Parameter, Reaction, Event, \
                      EventTrigger, EventAssignment, RateRule, \
                      AssignmentRule, FunctionDefinition, TimeSpan
from gillespy2.solvers import TauHybridCSolver, TauHybridSolver


class TestHybridEventRound(unittest.TestCase):

    def create_model(self, discrete=True):
        model = Model(name="test_issue710_round_event_assign")

        # Parameters
        k1 = Parameter(name="k1", expression="1e-05")
        model.add_parameter([k1])

        # Variables
        if discrete:
            A = Species(name="A", initial_value=5, mode='discrete')
        else:
            A = Species(name="A", initial_value=5, mode="continuous")
        model.add_species([A])

        # Reactions
        r1 = Reaction(name="r1", reactants={}, products={'A': 1}, rate="k1")
        model.add_reaction([r1])

        # Event Triggers
        e1_trig = EventTrigger(expression="t>10", initial_value=False, persistent=False)

        # Event Assignments
        e1_assign_1 = EventAssignment(variable="A", expression="A + 0.5")

        # Events
        model.add_event(Event(name="e1", trigger=e1_trig, assignments=[e1_assign_1], delay=None, priority="0", use_values_from_trigger_time=False))

        # Timespan
        tspan = TimeSpan(np.arange(0, 20.05, 0.05))
        model.timespan(tspan)
        return model

    def create_dynamic_model(self):
        model = Model(name="test_issue710_round_event_assign_switch")

        # Parameters
        k1 = Parameter(name="k1", expression="1e-05")
        model.add_parameter([k1])

        # Variables
        A = Species(name="A", initial_value=100, mode='dynamic', 
                    switch_min=50)
        model.add_species([A])

        # Reactions
        r1 = Reaction(name="r1", reactants={}, products={'A': 1}, rate="k1")
        model.add_reaction([r1])

        # Event Triggers
        t5_trig = EventTrigger(expression="t>5", initial_value=False, persistent=False)
        t10_trig = EventTrigger(expression="t>10", initial_value=False, persistent=False)
        t15_trig = EventTrigger(expression="t>15", initial_value=False, persistent=False)

        # Event Assignments
        eassign_1 = EventAssignment(variable="A", expression="A + 0.33")
        eassign_2 = EventAssignment(variable="A", expression="A - 90")

        # Events
        model.add_event(Event(name="e1", trigger=t5_trig, 
                             assignments=[eassign_1], 
                             delay=None, priority="0", 
                             use_values_from_trigger_time=False))
        model.add_event(Event(name="e2", trigger=t10_trig, 
                             assignments=[eassign_2], 
                             delay=None, priority="0", 
                             use_values_from_trigger_time=False))
        model.add_event(Event(name="e3", trigger=t15_trig, 
                             assignments=[eassign_1], 
                             delay=None, priority="0", 
                             use_values_from_trigger_time=False))

        # Timespan
        tspan = TimeSpan(np.arange(0, 20.05, 0.05))
        model.timespan(tspan)
        return model

    def setUp(self):
        self.d_model = self.create_model(discrete=True)
        self.c_model = self.create_model(discrete=False)

        self.dyn_model  = self.create_dynamic_model()

    def test_dynamic_event_rounding(self):
        expected_times_values = {
            99:  100.0,
            102: 100.33,
            202: 10.0,
            302: 10.0,
        }
        number_of_trajectories = 3

        for sname, sclass in {'TauHybridCSolver':TauHybridCSolver,'TauHybridSolver':TauHybridSolver}.items():
            with self.subTest(f"Checking event assignment rounding for dynamic species with {sname} solver."):
                solver = sclass(model=self.dyn_model)
                results = solver.run(number_of_trajectories=number_of_trajectories)
                for result in results:
                    for time_value, expected_value in expected_times_values.items():
                        self.assertAlmostEqual(result['A'][time_value], expected_value, places=2, msg=f"Simulation output at {time_value} is {result['A'][time_value]}, should be {expected_value}")


    def test_event_rounding(self):
        number_of_trajectories = 3
        for model, expected_last in {self.d_model: 6, self.c_model: 5.5}.items():
            expected_init = model.get_species("A").initial_value

            for sname, sclass in {'TauHybridCSolver':TauHybridCSolver,'TauHybridSolver':TauHybridSolver}.items():
                with self.subTest(f"Checking event assignment rounding for {model.listOfSpecies['A'].mode} species with {sname} solver."):
                    solver = sclass(model=model)
                    results = solver.run(number_of_trajectories=number_of_trajectories)
                    for result in results:
                        first_value = result["A"][0]
                        last_value =  result["A"][-1], 
                        self.assertAlmostEqual(first_value, expected_init, places=3, msg=f"Simulation output begins with {first_value}, should be {expected_init}")
                        self.assertAlmostEqual(result["A"][-1], expected_last, places=3, msg=f"Simulation output ends with {last_value}, should be {expected_last}")

