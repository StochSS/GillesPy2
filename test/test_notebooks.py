
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
import os
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor


class TestNotebooks(unittest.TestCase):
    
    def test_notebooks(self):
        FULL, MINIMAL = 0, 1
        test_set = MINIMAL
        notebooks = {}
        errors = {}
        ep = ExecutePreprocessor(timeout=600, kernel_name='python3', allow_errors=True)

        if test_set == FULL:
            test_dir = os.path.join('..', 'examples')
        elif test_set == MINIMAL:
            test_dir = os.path.join('..', os.path.join('examples', 'StartingModels'))

        for root, dirs, files in os.walk(test_dir):
            for file in files:
                if file.endswith(".ipynb"):
                     with open(os.path.join(root, file)) as f:
                        print('Reading {}...'.format(file))
                        notebooks[file] = nbformat.read(f, as_version=nbformat.NO_CONVERT)

        for file, nb in notebooks.items():
            with self.subTest(msg=file):
                try:
                    ep.preprocess(nb, {'metadata': {'path': root}})
                except Exception as err:
                    errors[file] = err
                self.assertFalse(bool(errors))


if __name__ == '__main__':
    unittest.main()

