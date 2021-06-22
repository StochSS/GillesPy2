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
            self.obj_dir = self.output_dir.joinpath("obj").resolve()
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

        self.output_file = "Simulation.out"
        if os.name == "nt":
            self.output_file = "Simulation.exe"

        self.output_file = Path(self.output_dir, self.output_file)

    def prebuild(self):
        self.__execute("prebuild")

    def build_simulation(self, simulation_name: str, **kwargs):
        self.__execute(simulation_name, **kwargs)

    def clean(self):
        self.__execute("clean")

    def clean_all(self):
        self.__execute("clean_all")

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
            f"- stderr: {result.stderr.decode('utf-8', errors='ignore')}\n")