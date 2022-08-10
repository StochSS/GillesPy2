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

from gillespy2.core import log
from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.cpp.ode_c_solver import ODECSolver
from gillespy2.solvers.cpp.tau_leaping_c_solver import TauLeapingCSolver
from gillespy2.solvers.cpp.tau_hybrid_c_solver import TauHybridCSolver

# Check to see if we're missing any dependencies.
from .build.build_engine import BuildEngine
missing_deps = BuildEngine.get_missing_dependencies()

if len(missing_deps) > 0:
    log.warning(
        f"Unable to use C++ optimized solvers due to one or more missing dependencies: {missing_deps}. "
        "The performance of this package can be significantly increased if you install/configure "
        "these on your machine."
    )

__all__ = ['SSACSolver', 'ODECSolver', 'TauLeapingCSolver', 'TauHybridCSolver']
