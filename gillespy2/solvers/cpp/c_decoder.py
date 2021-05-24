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
    def __init__(self):
        super(BasicSimDecoder, self).__init__()
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
        return super().get_output()

    def get_live_output(self):
        return super().get_live_output()

if __name__ == "__main__":
    test = BasicSimDecoder()
    print(test)
