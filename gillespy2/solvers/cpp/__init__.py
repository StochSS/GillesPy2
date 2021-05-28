from gillespy2.solvers.cpp import build
from gillespy2.core import logging
from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.cpp.ode_c_solver import ODECSolver
from gillespy2.solvers.cpp.tau_leaping_c_solver import TauLeapingCSolver
from gillespy2.solvers.cpp.variable_ssa_c_solver import VariableSSACSolver

# Check to see if we're missing any dependencies.
from .build.build_engine import BuildEngine
missing_deps = BuildEngine.get_missing_dependencies()

if len(missing_deps) > 0:
    logging.warn(
        f"Unable to use C++ optimized SSA due to one or more missing dependencies: {missing_deps}. "
        "The performance of this package can be significantly increased if you install/configure "
        "these on your machine."
    )

__all__ = ['SSACSolver', 'VariableSSACSolver']
__all__ = ['SSACSolver', 'ODECSolver', 'TauLeapingCSolver', 'VariableSSACSolver']
