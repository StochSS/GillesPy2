import unittest
import tempfile
import os
from gillespy2.core.model import import_SBML
from gillespy2 import ODESolver


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
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(response.read())
        tmp.close()
        sbml_model, errors = import_SBML(tmp.name)
        os.remove(tmp.name)
        sbml_results = sbml_model.run(solver=ODESolver)


if __name__ == '__main__':
    unittest.main()
