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
from collections.abc import Iterator

import numpy as np

from gillespy2.core.jsonify import Jsonify
from .gillespyError import TimespanError

class TimeSpan(Iterator, Jsonify):
    """
    Model timespan that describes the duration to run the simulation and at which timepoint to sample
    the species populations during the simulation.

    :param items: Evenly-spaced list of times at which to sample the species populations during the simulation.
            Best to use the form np.linspace(<start time>, <end time>, <number of time-points, inclusive>)
    :type items: list, tuple, range, or numpy.ndarray

    :raises TimespanError: items is an invalid type.
    """
    def __init__(self, items):
        if isinstance(items, np.ndarray):
            self.items = items
        elif isinstance(items, (list, tuple, range)):
            self.items = np.array(items)
        else:
            raise TimespanError("Timespan must be of type: list, tuple, range, or numpy.ndarray.")

        self.validate()

    def __str__(self):
        return self.items.__str__()

    def __eq__(self, o):
        return self.items.__eq__(o).all()

    def __getitem__(self, key):
        return self.items.__getitem__(key)

    def __iter__(self):
        return self.items.__iter__()

    def __len__(self):
        return self.items.__len__()

    def __next__(self):
        return self.items.__next__()

    @classmethod
    def linspace(cls, t=20, num_points=None):
        """
        Creates a timespan using the form np.linspace(0, <t>, <num_points, inclusive>).

        :param t: End time for the simulation.
        :type t: float | int

        :param num_points: Number of sample points for the species populations during the simulation.
        :type num_points: int

        :returns: Timespan for the model.
        :rtype: gillespy2.TimeSpan

        :raises TimespanError: t or num_points are None, <= 0, or invalid type.
        """
        if t is None or not isinstance(t, (int, float)) or t <= 0:
            raise TimespanError("t must be a positive float or int.")
        if num_points is not None and (not isinstance(num_points, int) or num_points <= 0):
            raise TimespanError("num_points must be a positive int.")

        if num_points is None:
            num_points = int(t / 0.05) + 1
        items = np.linspace(0, t, num_points)
        return cls(items)

    @classmethod
    def arange(cls, increment, t=20):
        """
        Creates a timespan using the form np.arange(0, <t, inclusive>, <increment>).

        :param increment: Distance between sample points for the species populations during the simulation.
        :type increment: float | int

        :param t: End time for the simulation.
        :type t: float | int

        :returns: Timespan for the model.
        :rtype: gillespy2.TimeSpan

        :raises TimespanError: t or increment are None, <= 0, or invalid type.
        """
        if t is None or not isinstance(t, (int, float)) or t <= 0:
            raise TimespanError("t must be a positive floar or int.")
        if not isinstance(increment, (float, int)) or increment <= 0:
            raise TimespanError("increment must be a positive float or int.")

        items = np.arange(0, t + increment, increment)
        return cls(items)

    def validate(self):
        """
        Validate the models time span

        :raises TimespanError: Timespan is an invalid type, empty, not uniform, contains a single \
                               repeated value, or contains a negative initial time.
        """
        if not isinstance(self.items, np.ndarray):
            if not isinstance(self.items, (list, tuple, range)):
                raise TimespanError("Timespan must be of type: list, tuple, range, or numpy.ndarray.")
            self.items = np.array(self.items)

        if len(self.items) == 0:
            raise TimespanError("Timespans must contain values.")
        if self.items[0] < 0:
            raise TimespanError("Simulation must run from t=0 to end time (t must always be positive).")

        first_diff = self.items[1] - self.items[0]
        other_diff = self.items[2:] - self.items[1:-1]
        isuniform = np.isclose(other_diff, first_diff).all()

        if not isuniform:
            raise TimespanError("StochKit only supports uniform timespans.")
        if first_diff == 0 or np.count_nonzero(other_diff) != len(other_diff):
            raise TimespanError("Timespan can't contain a single repeating value.")
