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
import numpy
import unittest

from gillespy2.core.timespan import TimeSpan
from gillespy2.core.gillespyError import TimespanError

class TestTimeSpan(unittest.TestCase):
    '''
    ################################################################################################
    Unit tests for gillespy2.TimeSpan.
    ################################################################################################
    '''
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

    def test_constructor__invalid_type(self):
        """ Test the TimeSpan constructor with an invalid data structure type. """
        test_tspan = set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
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

    def test_linspace__no_num_points(self):
        """ Test TimeSpan.linspace without passing num_points. """
        tspan = TimeSpan.linspace(t=30)
        self.assertEqual(tspan, numpy.linspace(0, 30, int(30 / 0.05) + 1))

    def test_linspace__no_args(self):
        """ Test TimeSpan.linspace without passing any args. """
        tspan = TimeSpan.linspace()
        self.assertEqual(tspan, numpy.linspace(0, 20, 401))

    def test_linspace__t_less_than_1(self):
        """ Test TimeSpan.linspace with t<1. """
        test_values = [0, -1, -2, -5, -10]
        for test_val in test_values:
            with self.subTest(t=test_val):
                with self.assertRaises(TimespanError):
                    tspan = TimeSpan.linspace(t=test_val, num_points=301)

    def test_linspace__num_points_less_than_1(self):
        """ Test TimeSpan.linspace with num_points<1. """
        test_values = [0, -1, -2, -5, -10]
        for test_val in test_values:
            with self.subTest(num_points=test_val):
                with self.assertRaises(TimespanError):
                    tspan = TimeSpan.linspace(t=30, num_points=test_val)

    def test_linspace__t_is_none(self):
        """ Test TimeSpan.linspace with t=None. """
        with self.assertRaises(TimespanError):
            tspan = TimeSpan.linspace(t=None, num_points=401)

    def test_linspace__invalid_t_type(self):
        """ Test TimeSpan.linspace with invalid t type. """
        with self.assertRaises(TimespanError):
            tspan = TimeSpan.linspace(t=[20.5], num_points=401)

    def test_linspace__invalid_num_points_type(self):
        """ Test TimeSpan.linspace with invalid num_points type. """
        with self.assertRaises(TimespanError):
            tspan = TimeSpan.linspace(t=20, num_points=40.1)

    def test_arange(self):
        """ Test TimeSpan.arange. """
        tspan = TimeSpan.arange(0.1, t=30)
        self.assertEqual(tspan, numpy.arange(0, 30.1, 0.1))

    def test_arange__no_t(self):
        """ Test TimeSpan.arange. """
        tspan = TimeSpan.arange(0.1)
        self.assertEqual(tspan, numpy.arange(0, 20.1, 0.1))

    def test_arange__t_less_than_1(self):
        """ Test TimeSpan.arange with t<1. """
        test_values = [0, -1, -2, -5, -10]
        for test_val in test_values:
            with self.subTest(t=test_val):
                with self.assertRaises(TimespanError):
                    tspan = TimeSpan.arange(0.1, t=test_val)

    def test_arange__num_points_less_than_1(self):
        """ Test TimeSpan.arange with increment<1. """
        test_values = [0, -1, -2, -5, -10]
        for test_val in test_values:
            with self.subTest(imcrement=test_val):
                with self.assertRaises(TimespanError):
                    tspan = TimeSpan.arange(test_val, t=30)

    def test_arange__t_is_none(self):
        """ Test TimeSpan.arange with t=None. """
        with self.assertRaises(TimespanError):
            tspan = TimeSpan.arange(0.1, t=None)

    def test_arange__invalid_t_type(self):
        """ Test TimeSpan.arange with invalid t type. """
        with self.assertRaises(TimespanError):
            tspan = TimeSpan.arange(0.05, t=[20.5])

    def test_arange__invalid_increment(self):
        """ Test TimeSpan.arange with invalid increment type. """
        with self.assertRaises(TimespanError):
            tspan = TimeSpan.arange("0.05", t=20)

    def test_validate__list(self):
        """ Test TimeSpan.validate with list data structure. """
        test_tspan = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        tspan = TimeSpan(test_tspan)
        tspan.items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        tspan.validate()
        self.assertEqual(tspan, numpy.array(test_tspan))

    def test_validate__tuple(self):
        """ Test TimeSpan.validate with tuple data structure. """
        test_tspan = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        tspan = TimeSpan(test_tspan)
        tspan.items = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        tspan.validate()
        self.assertEqual(tspan, numpy.array(test_tspan))

    def test_validate__range(self):
        """ Test TimeSpan.validate with range data structure. """
        test_tspan = range(11)
        tspan = TimeSpan(test_tspan)
        tspan.items = range(11)
        tspan.validate()
        self.assertEqual(tspan, numpy.array(test_tspan))

    def test_validate__invalid_type(self):
        """ Test TimeSpan.validate with an invalid data structure type. """
        test_tspan = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        tspan = TimeSpan(test_tspan)
        with self.assertRaises(TimespanError):
            tspan.items = set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
            tspan.validate()

    def test_validate__empty_timespan(self):
        """ Test TimeSpan.validate with an empty data structure. """
        test_tspan = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        tspan = TimeSpan(test_tspan)
        with self.assertRaises(TimespanError):
            tspan.items = []
            tspan.validate()

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
