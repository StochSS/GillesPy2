import io
import abc
from abc import ABC

class SimDecoder(ABC):
    """
    """
    def __init__(self):
        pass

    @abc.abstractmethod
    def read(self, output: io.BufferedReader):
        pass

    @abc.abstractmethod
    def get_output(self):
        pass

    @abc.abstractmethod
    def get_live_output(self):
        pass
