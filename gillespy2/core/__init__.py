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

import logging

from gillespy2.__version__ import __version__
from .assignmentrule import *
from .cleanup import *
from .events import *
from .functiondefinition import *
from .gillespyError import *
from .gillespySolver import *
from .model import *
from .parameter import *
from .raterule import *
from .reaction import *
from .results import *
from .sortableobject import *
from .species import *
from .timespan import TimeSpan

version = __version__

log = logging.getLogger("GillesPy2")
log.setLevel(logging.WARN)
log.propagate = False

if not log.handlers:
    _formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    _handler = logging.StreamHandler()
    _handler.setFormatter(_formatter)

    log.addHandler(_handler)

__all__ = [s for s in dir() if not s.startswith('_')]
