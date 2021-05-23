import shutil, tempfile
from pathlib import Path

from . import template_gen
from .make import Make
from gillespy2.core import Model

class BuildEngine():
    template_definitions_name = "template_definitions.h"

    def __init__(self, debug: bool = False):
        self.self_dir = Path(__file__).parent
        self.cpp_dir = self.self_dir.joinpath("../c_base").resolve()
        self.template_dir = self.cpp_dir.joinpath("template")
        self.debug = debug

        self.makefile = self.cpp_dir.joinpath("Makefile")


    def __enter__(self):
        if self.debug:
            self.temp_dir = Path(self.cpp_dir.joinpath("build"))

        else:
            self.temp_dir = Path(tempfile.mkdtemp())

        self.make = Make(self.makefile, self.temp_dir)

        return self

    def __exit__(self, type, value, traceback):
        """
        Delete the temp directory and all other associated build artifacts.
        """

        if self.debug:
            return

        shutil.rmtree(self.temp_dir)

    def prepare(self):
        """
        Copy the template directory to the temp directory.
        """

        # If object files haven't been compiled yet, go ahead and compile them with make.
        self.make.prebuild()

        shutil.copytree(self.template_dir, self.temp_dir, dirs_exist_ok=True)
        self.temp_dir.joinpath(self.template_definitions_name).unlink()

    def write_template(self, model: Model, variable=False):
        """
        Write model definitions to the template file. If the file does not exist, one is created.

        :param model: The model who's definitions will be written.
        :type model: gillespy2.Model

        :param variable: Set true to allow for non-constant parameter values.
        :type variable: bool
        """

        template_gen.write_template(self.temp_dir.joinpath(self.template_definitions_name), model, variable)

    def build_solver(self, solver_name: str):
        """
        Build the solver to the temp directory.

        :param solver_name: The name of the solver to build. For example, ODESolver.
        :type solver_name: str
        """

        self.make.build_solver(solver_name, template_dir=self.template_dir)

