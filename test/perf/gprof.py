import os
import subprocess
import time
from pathlib import Path

from gillespy2.core import Model
from gillespy2.solvers.cpp.c_solver import CSolver
from gillespy2.solvers.cpp.build.build_engine import BuildEngine

MAKEFILE_PATH = Path(os.path.dirname(__file__)).joinpath("Makefile")


def run_profiler(model: Model, solver: CSolver, trajectories=4, timesteps=101):
    """
    Run low-level performance tests on the specified solver.
    Profiling is performed with a system-level profiling tool.
    Only runtime performance is tested.

    Output is automatically parsed and formatted.
    """
    performance_result = 0

    # The simulation is built without running the solver directly.
    build: BuildEngine = BuildEngine()

    try:
        build.makefile = MAKEFILE_PATH
        build.prepare(model)

        # Prepare the location where performance data (gmon.out) is written.
        gmon_env = {
            "GMON_OUT_PREFIX": f"{str(build.output_dir)}/profile"
        }

        exe = build.build_simulation(simulation_name=solver.target)

        # Execute the simulation before running the profiler.
        process_args = [
            exe,
            "-trajectories", str(trajectories),
            "-timesteps", str(timesteps),
            "-end", "100"
        ]
        start = time.perf_counter()
        subprocess.check_call(args=process_args, stdout=subprocess.DEVNULL, env=gmon_env)
        stop = time.perf_counter()

        # Locate gprof profiling metadata in output directory.
        # Gprof will always append `.{pid}` to the end of {GMON_OUT_PREFIX}
        # As such, data is located by finding the file in the output dir
        #   which contains the matching prefix defined above.
        gprof_data: str = None
        for n in build.output_dir.iterdir():
            if n.stem == "profile":
                gprof_data = str(build.output_dir.joinpath(n.name))
                break
        if gprof_data is None:
            raise EnvironmentError("Profiler data was not found in the current environment; aborting")

        # Run gprof to process the profiler output.
        gprof_args = ["gprof", exe, gprof_data]
        perf_out = subprocess.check_output(gprof_args)
        # TODO: handle, parse, and process output from gprof
        print(perf_out)

        performance_result = round((stop - start) * 1000, 6)
    finally:
        build.clean()

    return performance_result
