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

import io
import abc
import numpy
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
    def read(self, output: io.BufferedReader):
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

    def read(self, output: io.BufferedReader):
        bytes_read = 0
        page_size = self.__read_next(output)
        while page_size > 0 and not output.closed:
            page_size = self.__read_next(output)
            bytes_read += page_size
        return bytes_read

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
