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

import unittest
import tempfile
import os
from example_models import RobustModel
from gillespy2.core.model import import_SBML, export_SBML
from gillespy2.core import gillespyError
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
        model = RobustModel()
        path = export_SBML(model, filename=tmp.name)
        document = libsbml.readSBML(path)
        if document.getNumErrors() > 0:
            error_id = document.getError(0).getErrorId()
            file_errors = [libsbml.XMLFileUnreadable, libsbml.XMLFileOperationError]
            assert error_id not in file_errors

    def test_sbml_model_without_propensities(self):
        try:
            import libsbml
        except ImportError:
            return

        model_path = os.path.join(os.path.dirname(__file__), "assets", "model_without_propensities.xml")

        try:
            import_SBML(model_path)
        except gillespyError.InvalidModelError:
            return

        self.fail()


if __name__ == '__main__':
    unittest.main()
