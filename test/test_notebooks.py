
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
import os
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor


class TestNotebooks(unittest.TestCase):
    
    def test_notebooks(self):
        errors = {}
        ep = ExecutePreprocessor(timeout=600, kernel_name='python3', allow_errors=True)
        for root, dirs, files in os.walk("../examples/"):
            for file in files:
                if file.endswith(".ipynb"):
                     with open(os.path.join(root, file)) as f:
                        print('Executing {}...'.format(file))
                        nb = nbformat.read(f, as_version=nbformat.NO_CONVERT)
                        try:
                            ep.preprocess(nb, {'metadata': {'path': root}})
                        except Exception as err:
                            print('Error executing the notebook "{}".\n\n'.format(file))
                            errors[file] = err
        for fname, err in errors.items():
            if len(err.__str__()) > 500:
                print('{}:\n{}\n...\n{}'.format(fname, err.__str__()[:251], err.__str__()[-251:]))
            else:
                print('{}:\n{}'.format(fname, err))

