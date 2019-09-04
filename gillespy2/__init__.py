"""
GillesPy2: Modeling toolkit for biochemical simulation
"""

import sys
if (sys.version_info < (3,0)):
    raise Exception("GillesPy2 only works in Python 3.0 and higher")

from .__version__ import __version__, __title__, __description__, __url__
from .__version__ import __author__, __email__
from .__version__ import __license__, __copyright__

from gillespy2.core import *
import gillespy2.sbml
import gillespy2.solvers
import gillespy2.example_models
