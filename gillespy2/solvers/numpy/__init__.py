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

import logging
_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_handler = logging.StreamHandler()
_handler.setFormatter(_formatter)
log = logging.getLogger()
log.setLevel(logging.WARN)
log.addHandler(_handler)

try:
    import numpy as np
    can_use_numpy = True
    from .ssa_solver import NumPySSASolver
    from .ode_solver import ODESolver
    from .tau_leaping_solver import TauLeapingSolver
    from .tau_hybrid_solver import TauHybridSolver
    log.debug("Successful Import of NumPy solvers.")

except Exception as e:
    log.warn(" Unable to use NumPy: {0}. The performance of this package can be significantly increased if you install NumPy.".format(e))
    can_use_numpy = False


__all__ = ['NumPySSASolver', 'ODESolver', 'TauLeapingSolver', 'TauHybridSolver'] if can_use_numpy else []
