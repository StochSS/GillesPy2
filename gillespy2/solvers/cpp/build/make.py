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
import subprocess

from pathlib import Path

from gillespy2.core import log
from gillespy2.core import gillespyError

class Make():
    def __init__(self, makefile: str, output_dir: str, obj_dir: str = None):
        self.makefile = Path(makefile).resolve()
        self.self_dir = Path(__file__).parent

        self.cbase_dir = self.self_dir.joinpath("../c_base").resolve()
        self.output_dir = Path(output_dir).resolve()

        # obj_dir is, presumably, only supplied if caching is enabled.
        # If not supplied, it should be assumed ethereal and cleaned up
        #  with the rest of the tmp directory.
        if obj_dir is None:
            self.obj_dir = self.output_dir.joinpath("gillespy2_obj").resolve()
        else:
            self.obj_dir = Path(obj_dir).resolve()

        for path in [self.output_dir, self.obj_dir]:
            if path.is_file():
                raise gillespyError.DirectoryError(
                    f"Error while attempting to open directory:\n"
                    f"- {path} is actually a file."
                )

            if not path.is_dir():
                path.mkdir()

        self.output_file = "GillesPy2_Simulation.out"
        if os.name == "nt":
            self.output_file = "GillesPy2_Simulation.exe"

        self.output_file = Path(self.output_dir, self.output_file)

    def prebuild(self):
        self.__execute("build")

    def build_simulation(self, simulation_name: str, **kwargs):
        self.__execute(simulation_name, **kwargs)

    def __execute(self, target: str, **kwargs):        
        # Default make arguments.
        args_dict = {
            "cbase_dir": str(self.cbase_dir.resolve()),
            "obj_dir": str(self.obj_dir.resolve()),
            "output_dir": str(self.output_dir.resolve()),
            "output_file": str(self.output_dir.joinpath(self.output_file).resolve())
        }

        # Overwrite keys supplied in **kwargs.
        for key, value in kwargs.items():
            args_dict[key] = value

        # Create the emake arguments. Note, arguments are UPPERCASE.
        make_args = [(f"{key.upper()}={value}") for key, value in args_dict.items()]

        # Create the make command.
        make_cmd = ["make", "-C", str(self.cbase_dir), "-f", str(self.makefile), target] + make_args

        try:
            result = subprocess.run(make_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        except KeyboardInterrupt:
            log.warning(f"Makefile was interrupted during execution of target: '{target}', unexpected behavior may occur.")

        if result.returncode == 0:
            return

        raise gillespyError.BuildError(f"Error encountered during execution of Makefile target: '{target}'.\n"
            f"Return code: {result.returncode}"
            f"- stdout: {result.stdout.decode('utf-8', errors='ignore')}\n"
            f"- stderr: {result.stderr.decode('utf-8', errors='ignore')}\n"
            f"- make_cmd: {make_cmd}\n"
            f"- os.listdir({os.path.join(self.cbase_dir,'template')}): {os.listdir(os.path.join(self.cbase_dir,'template'))}\n")
