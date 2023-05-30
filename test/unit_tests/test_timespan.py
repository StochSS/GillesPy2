# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2023 GillesPy2 developers.

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
''' Unit tests module for gillespy2.TimeSpan. '''
import time
import unittest
from datetime import datetime

import numpy

from gillespy2 import TimeSpan
from gillespy2 import TimespanError

class TestTimeSpan(unittest.TestCase):
    ''' Unit tests class for gillespy2.TimeSpan. '''
    def setUp(self):
        """ Setup a clean valid timespan for testing. """
        self.tspan = TimeSpan.linspace(t=100, num_points=101)

    def test_constructor(self):
        """ Test the TimeSpan constructor. """
        test_tspan = numpy.linspace(0, 20, 401)
        tspan = TimeSpan(test_tspan)
        self.assertEqual(tspan, test_tspan)

    def test_constructor__valid_data_structures(self):
        """ Test the TimeSpan constructor with valid data structures. """
        test_tspans = [
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            (1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
            range(11)
        ]
        for test_tspan in test_tspans:
            with self.subTest(tspan=test_tspan, tspan_type=type(test_tspan)):
                tspan = TimeSpan(test_tspan)
                self.assertEqual(tspan, numpy.array(test_tspan))

    def test_constructor__invalid_items(self):
        """ Test the TimeSpan constructor with an invalid data structure type. """
        test_tspans = [
            None, "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]", 20, 50.5,
            set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        ]
        for test_tspan in test_tspans:
            with self.subTest(items=test_tspan):
                with self.assertRaises(TimespanError):
                    TimeSpan(test_tspan)

    def test_linspace(self):
        """ Test TimeSpan.linspace. """
        tspan = TimeSpan.linspace(t=30, num_points=301)
        self.assertEqual(tspan, numpy.linspace(0, 30, 301))

    def test_linspace__no_t(self):
        """ Test TimeSpan.linspace without passing t. """
        tspan = TimeSpan.linspace(num_points=201)
        self.assertEqual(tspan, numpy.linspace(0, 20, 201))

    def test_linspace__invalid_t(self):
        """ Test TimeSpan.linspace with invalid t. """
        test_values = [None, "5", 0, -0.5, -1, -2, -5, -10, [20.5]]
        for test_val in test_values:
            with self.subTest(t=test_val):
                with self.assertRaises(TimespanError):
                    TimeSpan.linspace(t=test_val, num_points=301)

    def test_linspace__no_num_points(self):
        """ Test TimeSpan.linspace without passing num_points. """
        tspan = TimeSpan.linspace(t=30)
        self.assertEqual(tspan, numpy.linspace(0, 30, int(30 / 0.05) + 1))

    def test_linspace__invalid_num_points(self):
        """ Test TimeSpan.linspace with invalid num_points. """
        test_values = ["5", 0, -1, -2, -5, -10, 4.5, [40]]
        for test_val in test_values:
            with self.subTest(num_points=test_val):
                with self.assertRaises(TimespanError):
                    TimeSpan.linspace(t=30, num_points=test_val)

    def test_linspace__no_args(self):
        """ Test TimeSpan.linspace without passing any args. """
        tspan = TimeSpan.linspace()
        self.assertEqual(tspan, numpy.linspace(0, 20, 401))

    def test_arange(self):
        """ Test TimeSpan.arange. """
        tspan = TimeSpan.arange(0.1, t=30)
        self.assertEqual(tspan, numpy.arange(0, 30.1, 0.1))

    def test_arange__no_t(self):
        """ Test TimeSpan.arange. """
        tspan = TimeSpan.arange(0.1)
        self.assertEqual(tspan, numpy.arange(0, 20.1, 0.1))

    def test_arange__invalid_t(self):
        """ Test TimeSpan.arange with invalid t. """
        test_values = [None, "5", 0, -0.5, -1, -2, -5, -10, [20.5]]
        for test_val in test_values:
            with self.subTest(t=test_val):
                with self.assertRaises(TimespanError):
                    TimeSpan.arange(0.1, t=test_val)

    def test_arange__invalid_increment(self):
        """ Test TimeSpan.arange with invalid increment type. """
        test_values = [None, "0.05", 0, -1, -2, -5, -10, [0.05]]
        for test_val in test_values:
            with self.subTest(imcrement=test_val):
                with self.assertRaises(TimespanError):
                    TimeSpan.arange(test_val, t=30)

    def test_validate__list(self):
        """ Test TimeSpan.validate with list data structure. """
        test_tspan = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.tspan.items = test_tspan
        self.tspan.validate()
        self.assertIsInstance(self.tspan.items, numpy.ndarray)
        self.assertEqual(self.tspan, numpy.array(test_tspan))

    def test_validate__tuple(self):
        """ Test TimeSpan.validate with tuple data structure. """
        test_tspan = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        self.tspan.items = test_tspan
        self.tspan.validate()
        self.assertIsInstance(self.tspan.items, numpy.ndarray)
        self.assertEqual(self.tspan, numpy.array(test_tspan))

    def test_validate__range(self):
        """ Test TimeSpan.validate with range data structure. """
        test_tspan = range(11)
        self.tspan.items = test_tspan
        self.tspan.validate()
        self.assertIsInstance(self.tspan.items, numpy.ndarray)
        self.assertEqual(self.tspan, numpy.array(test_tspan))

    def test_validate__invalid_type(self):
        """ Test TimeSpan.validate with an invalid data structure type. """
        test_tspans = [
            None, "50", 20, 40.5, set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        ]
        for test_tspan in test_tspans:
            if test_tspan is not None:
                self.setUp()
            with self.subTest(items=test_tspan):
                with self.assertRaises(TimespanError):
                    self.tspan.items = test_tspan
                    self.tspan.validate()

    def test_validate__empty_timespan(self):
        """ Test TimeSpan.validate with an empty data structure. """
        test_tspans = [[], ()]
        for test_tspan in test_tspans:
            if test_tspan != []:
                self.setUp()
            with self.subTest(items=test_tspan):
                with self.assertRaises(TimespanError):
                    self.tspan.items = test_tspan
                    self.tspan.validate()

    def test_validate__all_same_values(self):
        """ Test TimeSpan.validate with an empty data structure. """
        test_tspan = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        tspan = TimeSpan(test_tspan)
        with self.assertRaises(TimespanError):
            tspan.items = [2, 2, 2, 2, 2, 2, 2, 2, 2]
            tspan.validate()

    def test_validate__negative_start(self):
        """ Test TimeSpan.validate with an initial time < 0. """
        test_tspan = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        tspan = TimeSpan(test_tspan)
        with self.assertRaises(TimespanError):
            tspan.items = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            tspan.validate()

    def test_validate__non_uniform_timespan(self):
        """ Test TimeSpan.validate with a non-uniform timespan. """
        test_tspan = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        tspan = TimeSpan(test_tspan)
        with self.assertRaises(TimespanError):
            tspan.items = [2, 1, 3, 4, 5, 6, 7, 8, 9, 10]
            tspan.validate()

    def test_comp_time_of_validate(self):
        """ Check the computation time of validate. """
        start = time.time()
        self.tspan.validate()
        tic = datetime.utcfromtimestamp(time.time() - start)
        msg = f"Total time to run validate on a timespan: {tic.strftime('%M mins %S secs %f msecs')}"
        print(f"\n<{'-'*88}>\n | {msg.ljust(84)} | \n<{'-'*88}>")
