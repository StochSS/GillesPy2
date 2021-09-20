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

from gillespy2.core import log
try:
    import pyximport
    import numpy as np
    pyximport.install(setup_args={'include_dirs': np.get_include()})
    from gillespy2.solvers.cython.cython_ssa_solver import CythonSSASolver
    can_use_cython = True
    log.debug("Successful Import of Cython solvers.")
except Exception as e:
    can_use_cython = False

__all__ = ['CythonSSASolver'] if can_use_cython else []
