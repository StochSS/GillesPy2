import unittest
import tempfile
import os
import ssl
from example_models import Example
from gillespy2.core.model import import_SBML, export_SBML
from gillespy2 import ODESolver


class TestSBML(unittest.TestCase):

    def test_sbml_conversion(self):
        try:
            import libsbml
        except ImportError:
            return

        try:
            from urllib2 import urlopen, URLError

        except ImportError:
            from urllib.request import urlopen
            from urllib.error import URLError

        sbml_file = 'https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000028.2?filename=BIOMD0000000028_url.xml'
        try:
            response = urlopen(sbml_file)
        except URLError:
            from certifi import where
            ctx = ssl.create_default_context()
            ctx.load_default_certs()
            ctx.load_verify_locations(where())
            response = urlopen(sbml_file, context=ctx)
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(response.read())
        tmp.close()
        sbml_model, errors = import_SBML(tmp.name)
        os.remove(tmp.name)
        sbml_results = sbml_model.run(solver=ODESolver)

    def test_sbml_export_conversion(self):
        try:
            import libsbml
        except ImportError:
            return

        tmp = tempfile.NamedTemporaryFile(delete=False)
        model = Example()
        path = export_SBML(model, filename=tmp.name)
        document = libsbml.readSBML(path)
        if document.getNumErrors() > 0:
            error_id = document.getError(0).getErrorId()
            file_errors = [libsbml.XMLFileUnreadable, libsbml.XMLFileOperationError]
            assert error_id not in file_errors


if __name__ == '__main__':
    unittest.main()
