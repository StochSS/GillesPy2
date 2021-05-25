import io
import abc
import numpy
from abc import ABC

class SimDecoder(ABC):
    """
    Abstract simulation decoder class.
    For solvers to handle output in a custom way,
      this decoder class will be implemented.
    Expects the output to be read from a buffered IO reader.
    """
    
    def __init__(self, trajectories: numpy.ndarray):
        """
        Constructor for simulation decoder class.
        Reads the output of an external simulation and generates a numpy array
          from the results.

        :param trajectories: 3D array to output simulation data to.
        :type trajectories: numpy.array
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

    @abc.abstractmethod
    def read(self, output: io.BufferedReader):
        pass

    @abc.abstractmethod
    def get_output(self):
        pass

    @abc.abstractmethod
    def get_live_output(self):
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
        """
        line = output.read().decode("utf-8")
        ln = len(line)
        if ln > 0:
            self.buffer.append(line)
        return ln

    def read(self, output: io.BufferedReader):
        """
        """
        bytes_read = 0
        page_size = self.__read_next(output)
        while page_size > 0 and not output.closed:
            page_size = self.__read_next(output)
            bytes_read += page_size
        return bytes_read

    def get_output(self):
        """
        Returns the fully-populated NumPy array containing the completed simulation data.
        Assumes that the subprocess has already completed.

        :return: Tuple containing the 3D NumPy array of results, and the time stopped.
        """
        stdout = "".join(self.buffer).split(",")
        # The last number written to stdout from C++ sim is always the stop time.
        time_stopped = stdout.pop()

        # Assumed layout of NumPy array:
        #  1D: index to each simulation trajectory
        #  2D: index to each timestep of that directory
        #  3D: index to each species of that timestep
        # Buffer is a flat 1D list, which gets mapped into the NumPy array.
        for traj_number, trajectory in enumerate(self.trajectories):
            # traj_i is this trajectory's offset into the 1D list.
            # Each trajectory has a "stride" equal to the total number of timesteps.
            traj_i = traj_number * self.num_timesteps

            for ts_number, timestep in enumerate(trajectory):
                # time_i is this timestep's offset into the 1D list.
                # Each timestep has a "stride" equal to the total number of timesteps.
                time_i = ts_number * self.num_species

                for spec_i in range(timestep.size):
                    # Output is a 1-dimensional list with the assumed layout.
                    # This species's offset is relative to the offset of the current timestep.
                    current_index = traj_i + time_i + spec_i
                    timestep[spec_i] = stdout[current_index]

        return self.trajectories, time_stopped

    def get_live_output(self):
        """
        Returns an iterator to iterate over the solver results as they come in.
        NOT YET IMPLEMENTED!
        """
        return super().get_live_output()
