import unittest
import tempfile
import os
from gillespy2.core import Model, import_SBML
from gillespy2.sbml import SBMLimport
from gillespy2.core.gillespyError import *
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
import numpy as np

class TestSBML(unittest.TestCase):

    def test_sbml_conversion(self):
        try:
            import libsbml
        except ImportError:
            return

        try:
            from urllib2 import urlopen

        except ImportError:
            from urllib.request import urlopen

        sbml_file = 'http://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD0000000028'
        response = urlopen(sbml_file)
        tmp = tempfile.NamedTemporaryFile(delete = False)
        tmp.write(response.read())
        tmp.close()
        sbml_model, errors = import_SBML(tmp.name)
        os.remove(tmp.name)
        sbml_results = sbml_model.run(solver=BasicODESolver)

if __name__ == '__main__':
    unittest.main()
