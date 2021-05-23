import logging
import os, subprocess

from pathlib import Path

from numpy import std
from gillespy2.core import gillespyError, logging

# cmd = ["make", "-C", self.output_directory, '-f', MAKE_FILE,
#         'ODESimulation', 'GILLESPY_C_ODE_DIR='+GILLESPY_C_ODE_DIR, 'CBASE_DIR='+CBASE_DIR,
#         'SUNDIALS_DIR='+SUNDIALS_DIR]

class Make():
    def __init__(self, makefile: str, output_dir: str):
        self.makefile = Path(makefile).resolve()
        self.self_dir = Path(__file__).parent

        self.cbase_dir = self.self_dir.joinpath("../c_base").resolve()
        self.obj_dir = self.cbase_dir.joinpath("obj").resolve()
        self.output_dir = Path(output_dir).resolve()

        if not self.output_dir.is_dir():
            self.output_dir.mkdir()

        if not self.obj_dir.is_dir():
            self.obj_dir.mkdir()

    def prebuild(self):
        self.__execute("prebuild")

    def build_solver(self, solver_name: str, **kwargs):
        self.__execute(solver_name, **kwargs)

    def clean(self):
        self.__execute("clean")

    def clean_all(self):
        self.__execute("clean_all")

    def __execute(self, target: str, **kwargs):        
        # Default make arguments.
        args_dict = {
            "cbase_dir": self.cbase_dir.resolve(),
            "obj_dir": self.obj_dir.resolve(),
            "output_dir": self.output_dir.resolve()
        }

        # Overwrite keys supplied in **kwargs.
        for key, value in kwargs.items():
            args_dict[key] = value

        # Create the emake arguments. Note, arguments are UPPERCASE.
        make_args = [(f"{key.upper()}={value}") for key, value in args_dict.items()]

        # Create the make command.
        make_cmd = ["make", "-C", self.cbase_dir, "-f", self.makefile, target] + make_args

        print(make_cmd)

        try:
            result = subprocess.run(make_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE)

        except KeyboardInterrupt:
            logging.warn(f"Makefile was interrupted during execution of target: '{target}', unexpected behavior may occur.")

        if result.returncode == 0:
            return

        raise gillespyError.BuildError(f"Error encountered during execution of Makefile target: '{target}'.\n"
            f"Return code: {result.returncode}"
            f"- stdout: {result.stdout.decode('utf-8')}\n")