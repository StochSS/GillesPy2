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

import os
import shutil
import tempfile

def cleanup_tempfiles():
    '''
    Cleanup all tempfiles in gillespy2 build.
    '''
    tempdir = tempfile.gettempdir()
    for file_obj in os.listdir(tempdir):
        if file_obj.startswith("gillespy2_build"):
            shutil.rmtree(os.path.join(tempdir, file_obj))
