import unittest
import gillespy2
from gillespy2 import TauHybridCSolver
import numpy as np

class EventFeatures(unittest.TestCase):
    class BaseEventModel(gillespy2.Model):
        def __init__(self, s1, s2, rate):
            super().__init__(name="BasicEventModel")

            self.s1 = gillespy2.Species(name="S1", initial_value=s1, mode="continuous")
            self.s2 = gillespy2.Species(name="S2", initial_value=s2, mode="continuous")
            self.rate = gillespy2.Parameter(name="k1", expression=rate)
            self.add_species([self.s1, self.s2])
            self.add_parameter([self.rate])
            self.add_reaction([
                gillespy2.Reaction(name="r1", rate=self.rate,
                                   reactants={self.s1: 1},
                                   products={self.s2: 1})
            ])

    def test_event_with_time_trigger(self):
        model = EventFeatures.BaseEventModel(s1=0, s2=0, rate=0.0)
        event = gillespy2.Event(name="ev1", assignments=[
            gillespy2.EventAssignment(variable=model.s1, expression="100.0"),
            gillespy2.EventAssignment(variable=model.rate, expression="1.0")
        ], trigger=gillespy2.EventTrigger(expression="t>5"))
        model.add_event(event)

        solver = TauHybridCSolver(model=model)
        result = model.run(solver=solver)[0]
        s1, s2 = result["S1"][-1], result["S2"][-1]

        self.assertGreater(s2, s1, "Expected S2 > S1")
        self.assertGreater(s1, 0.0, "Expected S1 > 0")
        self.assertAlmostEqual(s1 + s2, 100.0, places=1)

    def test_event_with_species_trigger(self):
        model = EventFeatures.BaseEventModel(s1=100, s2=0, rate=10.0)
        event = gillespy2.Event(name="ev1", assignments=[
            gillespy2.EventAssignment(variable=model.s1, expression="100.0"),
            gillespy2.EventAssignment(variable=model.rate, expression="0.0")
        ], trigger=gillespy2.EventTrigger(expression="S1<90"))
        model.add_event(event)

        solver = TauHybridCSolver(model=model)
        result = model.run(solver=solver)[0]
        s1, s2 = result["S1"][-1], result["S2"][-1]

        self.assertEqual(s1, 100, "Expected S1 == 100 (trigger set S1 to 100 and rate to 0")
        self.assertGreater(s2, 0, "Expected S2 > 0")
        self.assertFalse(np.any(result["S1"] <= 90.0), "Expected S1 > 90 for entire simulation")

    def test_delay_trigger_persistent(self):
        model = EventFeatures.BaseEventModel(s1=100, s2=0, rate=1.0)
        event1 = gillespy2.Event(name="ev1", assignments=[
            gillespy2.EventAssignment(variable=model.s1, expression="0"),
            gillespy2.EventAssignment(variable=model.s2, expression="0"),
            gillespy2.EventAssignment(variable=model.rate, expression="0.0")
        ], trigger=gillespy2.EventTrigger(expression="S1<60 and S2<S1", persistent=False), delay="t+1.0")
        event2 = gillespy2.Event(name="ev2", assignments=[
            gillespy2.EventAssignment(variable=model.s1, expression="200"),
            gillespy2.EventAssignment(variable=model.s2, expression="200"),
            gillespy2.EventAssignment(variable=model.rate, expression="0.0"),
        ], trigger=gillespy2.EventTrigger(expression="S2>90 and t<3.5", persistent=True), delay="1.0")
        model.add_event([event1, event2])

        solver = TauHybridCSolver(model=model)
        result = model.run(solver=solver)[0]
        s1, s2 = result["S1"][-1], result["S2"][-1]

        # If delay is working correctly:
        # * event1 is never triggered. event1 sets everything to 0.
        #   If event1 fires, event2 can never fire.
        # * event2 is triggered, setting everything to 100 (and rate to 0).
        self.assertNotIn(0, [s1, s2], "Non-persistent event fired unexpectedly")
        self.assertEqual(s1, 200, "Persistent event failed to fire")
        self.assertEqual(s2, 200, "Persistent event failed to fire")

    def test_trigger_priorities(self):
        model = EventFeatures.BaseEventModel(s1=100, s2=0, rate=1.0)
        event1 = gillespy2.Event(name="ev1", assignments=[
            gillespy2.EventAssignment(variable=model.s1, expression="100"),
            gillespy2.EventAssignment(variable=model.s2, expression="100"),
        ], trigger=gillespy2.EventTrigger(expression="S1 < 50"), priority="2*t*S1")
        event2 = gillespy2.Event(name="ev2", assignments=[
            gillespy2.EventAssignment(variable=model.s1, expression="0"),
        ], trigger=gillespy2.EventTrigger(expression="S1 < 50"), priority="t*S1")
        model.add_event([event1, event2])

        solver = TauHybridCSolver(model=model)
        result = model.run(solver=solver)[0]
        s1, s2 = result["S1"][-1], result["S2"][-1]

        # If priority is working correctly, event2 should ALWAYS fire before event1.
        # Proper result is S1 = 0, S2 = 100, so no further reactions are possible.
        self.assertEqual(s1, 0,   "Events fired in an incorrect order")
        self.assertEqual(s2, 100, "Events fired in an incorrect order")

    def test_use_values_from_trigger_time(self):
        model = EventFeatures.BaseEventModel(s1=100, s2=0, rate=1.0)
        event = gillespy2.Event(name="ev1", assignments=[
            gillespy2.EventAssignment(variable=model.s1, expression="S2"),
            gillespy2.EventAssignment(variable=model.rate, expression="0.0"),
        ], trigger=gillespy2.EventTrigger(expression="S1 < 60"), delay="1.5", use_values_from_trigger_time=True)
        model.add_event(event)

        solver = TauHybridCSolver(model=model)
        result = model.run(solver=solver)[0]
        s1, s2 = result["S1"][-1], result["S2"][-1]

        self.assertGreater(s2, s1, "Event assignment did not assign values from trigger time")

    def test_initial_values(self):
        model = EventFeatures.BaseEventModel(s1=0, s2=100.0, rate=1.0)

        event = gillespy2.Event(name="ev1", assignments=[
            gillespy2.EventAssignment(variable=model.s1, expression="S2/2"),
        ], trigger=gillespy2.EventTrigger(expression="S1==0", initial_value=False))
        model.add_event(event)

        solver = TauHybridCSolver(model=model)
        result = model.run(solver=solver)

        s1, s2 = result["S1"][-1], result["S2"][-1]
        self.assertAlmostEqual(s1 + s2, 150.0, places=1, msg="Event assignment assigned incorrect value")
        self.assertGreater(s2, 100, "Event with initial condition did not fire")
        self.assertEqual(result["S1"][0], 50, "Event assignment with initial condition failed to fire at t=0")
