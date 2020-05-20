import logging

from gillespy2.core import gillespyError
from gillespy2.core.Model import *
from gillespy2.core.Species import Species
from gillespy2.core.Reaction import Reaction
from gillespy2.core.Parameter import Parameter
from gillespy2.core.RateRule import RateRule
from gillespy2.core.Utilities import *
from gillespy2.core.FunctionDefinition import FunctionDefinition
from gillespy2.core.AssignmentRule import AssignmentRule
from gillespy2.core.SortableObject import SortableObject
from gillespy2.core.gillespySolver import GillesPySolver
from gillespy2.core.events import *
from gillespy2.__version__ import __version__


_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_handler = logging.StreamHandler()
_handler.setFormatter(_formatter)

version = __version__
log = logging.getLogger()
log.setLevel(logging.WARN)
log.addHandler(_handler)

__all__ = [s for s in dir() if not s.startswith('_')]
