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

import unittest
import tempfile
import os
import re
from example_models import create_robust_model
from gillespy2.core.model import import_SBML, export_SBML
from gillespy2.core import gillespyError
from gillespy2 import ODESolver


class TestSBML(unittest.TestCase):

    converted_symbols = {
        "ln": "log",
    }
    test_formula_name = "ar_t_ave"

    def test_sbml_conversion(self):
        try:
            import libsbml
        except ImportError:
            return

        model_path = os.path.join(os.path.dirname(__file__), "assets", "test_sbml.xml")
        sbml_model, errors = import_SBML(model_path)
        solver = ODESolver(model=sbml_model)
        sbml_results = sbml_model.run(solver=solver, increment=0.05, t=20)

    def test_sbml_string_replace(self):
        """
        Ensure that non-supported symbols are converted to supported symbols (such as ln -> log)
        """
        try:
            import libsbml

            # Model is modified from BIOMD0000000012, with species names modified to force possible collisions.
            # For example, species name `PY` -> `xlnx` which, if not implemented correctly,
            #   would cause a collision with `ln` and result in the name `xlogx`, which is wrong.
            model_path = os.path.join(os.path.dirname(__file__), "assets", "math_test_model.xml")
            sbml_model, errors = import_SBML(model_path)
            test_rule = sbml_model.listOfAssignmentRules.get(TestSBML.test_formula_name).formula

            for before, after in TestSBML.converted_symbols.items():
                with self.subTest("Testing SBML conversion", before=before, after=after):
                    self.assertTrue(f"x{before}x" in sbml_model.listOfSpecies,
                                    "Species name 'xlnx' was unexpectedly modified during SBML import")
                    self.assertIsNone(re.search(rf'\b{before}\b', test_rule),
                                      f"Expected SBML formula '{test_rule}' to remove all instances of '{before}'")
                    self.assertIsNotNone(re.search(rf'\b{after}\b', test_rule),
                                         f"Expected SBML formula '{test_rule} to convert all instances of '{after}'")
        except ImportError:
            return

    def test_sbml_export_conversion(self):
        try:
            import libsbml
        except ImportError:
            return

        tmp = tempfile.NamedTemporaryFile(delete=False)
        model = create_robust_model()
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
        except gillespyError.ModelError:
            return

        self.fail()


if __name__ == '__main__':
    unittest.main()
