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

import io
import abc
import numpy
from collections import deque
from typing import Callable
from abc import ABC

class SimDecoder(ABC):
    """
    Abstract simulation decoder class.
    For solvers to handle output in a custom way, this decoder class will be implemented.
    Expects the output to be read from a buffered IO reader.

    :param trajectories: 3D array to output simulation data to.
    :type trajectories: numpy.array
    """
    
    def __init__(self, trajectories: numpy.ndarray):
        """
        Constructor for simulation decoder class.
        Reads the output of an external simulation and generates a numpy array from the results.
        """
        # Make sure that it's actually a numpy array
        if not isinstance(trajectories, numpy.ndarray):
            raise io.UnsupportedOperation("Could not construct decoder: only NumPy arrays are supported")
        # ...and also that it's actually 3-dimensional
        if not len(trajectories.shape) == 3:
            raise io.UnsupportedOperation("Could not construct decoder: supplied NumPy array must be 3-dimensional")

        # Provided trajectory assumed valid
        self.trajectories = trajectories
        self.num_trajectories, self.num_timesteps, self.num_species = trajectories.shape
        self.bytes_read = 0

    @classmethod
    def create_default(cls, num_trajectories: int, num_timesteps: int, num_species: int) -> "SimDecoder":
        """
        Creates a new instance of the calling class, using a NumPy array with a predefined shape.
        Calling this method is preferred over calling the constructor directly.

        :param num_trajectories: Number of trajectories expected in the simulation output.
        :type num_trajectories: int

        :param num_timesteps: Number of timesteps expected in the simulation output.
        :type num_timesteps: int

        :param num_species: Number of species expected in the simulation output.
        :type num_species: int

        :returns: An instance of the decoder object, automatically populated with a valid output array.
        """
        return cls(numpy.zeros((num_trajectories, num_timesteps, num_species + 1)))

    @abc.abstractmethod
    def read(self, output: io.BufferedReader, **kwargs):
        """
        Accepts a buffered reader from stdout of a subprocess.
        Contents of the given reader are processed and made available through get_output().

        Blocks until the output of the buffered reader has been read completely.

        :param output: Reader provided from the stdout member of an open Popen class.
        :type output: io.BufferedReader
        """
        pass

    @abc.abstractmethod
    def get_output(self) -> "tuple[numpy.ndarray, int]":
        """
        Returns the fully-populated NumPy array containing the completed simulation data.
        Assumes that the subprocess has already completed.

        :returns: Tuple containing the 3D NumPy array of results, and the time stopped.
        """
        pass

class BasicSimDecoder(SimDecoder):
    """
    Simple decoder which returns the results as a complete string.
    """
    def __init__(self, trajectories: numpy.ndarray):
        super(BasicSimDecoder, self).__init__(trajectories)
        self.buffer = []

    def __read_next(self, output: io.BufferedReader):
        """
        Reads the next block from the simulation output.
        Returns the length of the string read.

        :param output: Reader provided from the stdout member of an open Popen object.
        :type output: io.BufferedReader

        :returns: Integer representing the number of characters read.
        """
        line = output.read().decode("utf-8")
        ln = len(line)
        if ln > 0:
            self.buffer.append(line)
        return ln

    def read(self, output: io.BufferedReader, **kwargs):
        page_size = self.__read_next(output)
        self.bytes_read = page_size
        while page_size > 0 and not output.closed:
            page_size = self.__read_next(output)
            self.bytes_read += page_size
        return self.bytes_read

    def get_output(self):
        stdout = "".join(self.buffer).split(",")
        # The last number written to stdout from C++ sim is always the stop time.
        time_stopped = stdout.pop()
        time_stopped = int(time_stopped) if time_stopped.isdigit() else 0

        # Assumed layout of NumPy array:
        #  1D: index to each simulation trajectory
        #  2D: index to each timestep of that directory
        #  3D: index to each species of that timestep
        # Buffer is a flat 1D list, which gets mapped into the NumPy array.
        for entry_i, entry in enumerate(stdout):
            spec_num = entry_i % self.num_species

            # Each timestep has a "stride" equal to the total number of species.
            entry_i //= self.num_species
            ts_num = entry_i % self.num_timesteps

            # Each trajectory has a "stride" equal to the total number of timesteps.
            traj_num = entry_i // self.num_timesteps

            self.trajectories[traj_num][ts_num][spec_num] = entry

        return self.trajectories, time_stopped


class IterativeSimDecoder(SimDecoder):
    """
    Output decoder for processing the output at regular intervals.
    An IterativeSimDecoder object can be provided a callback.
    Each time a new timestep is received and processed, the callback is invoked,
      and is provided with the timestep value and output time.

    Intended for handling just-in-time output.
    """

    def __init__(self, trajectories: numpy.ndarray, callback=None):
        super(IterativeSimDecoder, self).__init__(trajectories)
        self.callback = callback if callback is not None else lambda *args: None
        self.end_time = 0

    def with_callback(self, callback: "Callable[[int, numpy.ndarray], None]") -> "IterativeSimDecoder":
        """
        Provide a callback function to be invoked on each timestep.
        Accepts a function with the signature:

        def callback(float, numpy.ndarray)

        The first value (float) is the time value for the given timestep entry.
        The second value (NumPy array) is the simulation state for that timestep.
        Return values are ignored.

        :param callback: Function to be invoked on each timestep.
        :type callback: Callable[[int, numpy.ndarray], None]

        :returns: Pass-through for IterativeSimDecoder.
        """
        if not callable(callback):
            raise ValueError(f"Expected function as callback: got {type(callback)}")
        self.callback = callback
        return self

    def read(self, output: io.BufferedReader, page_size=256, **kwargs):
        """
        Read and process output from the provided buffer, one timestep at a time.
        Any registered callbacks will be invoked on each iteration of output processing.

        Blocks until the output of the buffered reader has been read completely.

        :param output: Reader provided from the stdout member of an open Popen class.
        :type output: io.BufferedReader

        :param page_size: Suggested value for number of bytes to read from the simulation on each pass.
        Smaller values may result in more consistent callback times, at the cost of performance.
        Larger values may result in better overall performance, at the cost of sporadic callback times.
        :type page_size: int

        :returns: Total number of bytes read
        """
        if page_size < 1:
            page_size = 256
        traj_id, t = 0, 0
        current_timestep = deque()
        num_trajectories, num_timesteps, entries_per_timestep = self.trajectories.shape

        line = output.read(page_size).decode("ascii")
        carry_value = ""
        bytes_read = len(line)
        self.bytes_read = bytes_read
        while bytes_read > 0:
            # carry_value is carried over from the previous read.
            # It must be appended to the beginning of the next read.
            if carry_value != "":
                line = carry_value + line
            entries = line.split(",")

            # The last "value" is not guaranteed to be a complete number.
            # It may be the start of an incomplete value. Example:
            # [ current read | next read ]
            # [...12,13,14,1 | 5,16,17...]
            # The value at the "boundary" should be 15, but if we do not carry,
            #   then it will be misread as two separate values, 1 and 5.
            carry_value = entries.pop()

            # Process incoming text, one timestep entry at a time.
            current_timestep.extend(entries)
            while len(current_timestep) >= entries_per_timestep:
                self.trajectories[traj_id][t] = [current_timestep.popleft() for _ in range(entries_per_timestep)]
                # First value of each timestep is the current time of the simulation.
                # The remaining entries of each timestep are the output state at that time.
                self.callback(self.trajectories[traj_id][t][0],
                              self.trajectories[traj_id][t][1:])
                t += 1
                if t >= num_timesteps:
                    traj_id += 1
                    t = 0

            # Process the next output block.
            # An empty output block indicates that the simulation has ended.
            line = output.read(page_size).decode("ascii")
            bytes_read = len(line)
            self.bytes_read += bytes_read

        if carry_value != "":
            self.end_time = int(carry_value)
        elif len(current_timestep) > 0:
            self.end_time = int(current_timestep.popleft())
        else:
            self.end_time = 0

        return self.bytes_read

    def get_output(self) -> "tuple[numpy.ndarray, int]":
        # TODO: block get_output() call if waiting on read() to complete
        return self.trajectories, self.end_time
