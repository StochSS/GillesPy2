import logging
import os, subprocess

from gillespy2.core import gillespyError, logging

# cmd = ["make", "-C", self.output_directory, '-f', MAKE_FILE,
#         'ODESimulation', 'GILLESPY_C_ODE_DIR='+GILLESPY_C_ODE_DIR, 'CBASE_DIR='+CBASE_DIR,
#         'SUNDIALS_DIR='+SUNDIALS_DIR]

class Make():
    def __init__(self, makefile: str, output_dir: str):
        self.makefile = os.path.abspath(makefile)
        self.self_dir = os.path.dirname(os.path.abspath(__file__))

        self.cbase_dir = os.path.join(self.self_dir, "../c_base")
        self.obj_dir = os.path.join(self.cbase_dir, "obj")
        self.cache_dir = os.path.join(self.cbase_dir, "cache")
        self.output_dir = os.path.abspath(output_dir)

        if not os.path.isdir(self.cache_dir):
            os.makedirs(self.cache_dir)

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

    def prebuild(self):
        self.__execute("prebuild", obj_dir=self.cache_dir)

    def __execute(self, target: str, **kwargs):        
        # This is the default make arguments. This needs to be merged with kwargs, keys to uppercase, and then formatted.
        make_args = {
            "cbase_dir": self.cbase_dir,
            "obj_dir": self.obj_dir,
            "output_dir": self.output_dir
        }
        make_args = [(f"{key.upper()}={value}") for key, value in make_args.items()]

        # Create the make command.
        make_cmd = ["make", "-C", self.cbase_dir, "-f", self.makefile, target] + make_args

        try:
            result = subprocess.run(make_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE)

        except KeyboardInterrupt:
            logging.warn(f"Makefile was interrupted during execution of target: '{target}', unexpected behavior may occur.")

        if result.returncode == 0:
            return

        raise gillespyError.BuildError(f"Error encountered during execution of Makefile target: '{target}'.\n"
            f"Return code: {result.returncode}"
            f"- stdout: {result.stdout.decode('utf-8')}\n"
            f"- stderr: {result.stderr.decode('utf-8')}\n")

make = Make("gillespy2/solvers/cpp/c_base/Makefile", "gillespy2/solvers/cpp/c_base/output")
make.prebuild()