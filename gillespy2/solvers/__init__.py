from gillespy2.core.gillespyError import ExecutionError
from .cpp import *
from .numpy import *

__all__ = cpp.__all__ + numpy.__all__

if not __all__:
    raise ExecutionError("Your computer does not contain the minimum require"
                         "ments for running simulations using GillesPy2."
                         " Please install NumPy, or configure the C++ 'g++' "
                         "compiler on your machine.")
