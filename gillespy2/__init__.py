import sys
if (sys.version_info < (3,0)):
    raise Exception("GillesPy2 only works in Python 3.0 and higher")


from gillespy2.core import *
import gillespy2.sbml
import gillespy2.solvers
import gillespy2.example_models
