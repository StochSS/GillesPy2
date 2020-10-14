from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.cpp.variable_ssa_c_solver import VariableSSACSolver
from gillespy2.core import log

from gillespy2.solvers.utilities.cpp_support_test import cpp_support
can_use_cpp = cpp_support

__all__ = ['SSACSolver', 'VariableSSACSolver']
