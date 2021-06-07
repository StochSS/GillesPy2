import unittest
import tempfile
import os
from example_models import Example
from gillespy2.core.model import import_SBML, export_SBML
from gillespy2 import ODESolver


class TestSBML(unittest.TestCase):

    def test_sbml_conversion(self):
        try:
            import libsbml
        except ImportError:
            return

        model_path = os.path.join(os.path.dirname(__file__), "assets", "test_sbml.xml")
        sbml_model, errors = import_SBML(model_path)
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
