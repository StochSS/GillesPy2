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

import os
import re
import subprocess
import time
from pathlib import Path

from .performance_data import PerformanceData
from .performance_data import PerformanceEntry

from gillespy2.core import Model
from gillespy2.solvers.cpp.c_solver import CSolver
from gillespy2.solvers.cpp.build.build_engine import BuildLevel
from gillespy2.solvers.cpp.build.build_engine import BuildEngine


def parse_gprof_output(output: str):
    # For each line, attempt to parse it as 3-6 floats and a name.
    # This is how each performance entry is listed.
    def parse_line(block: str):
        block = re.split("\\s+", block.strip())

        # Check to make sure we have an expected number of entries.
        if len(block) < 4:
            return "", None, 0

        # Check to make sure each entry is parsable as a float
        elif not all(n.replace(".", "", 1).isnumeric() for n in block[:3]):
            return "", None, 0

        # 1: % spent on current function call
        # 2: cumulative time (time spent on this call + the above)
        percent, cumulative_t, t = block[:3]
        # Remove fully-numeric entries
        block = [b for b in block[3:] if not b.replace(".", "", 1).isnumeric()]
        call_name = "".join(block)

        result = PerformanceEntry(t=t, percent=percent)
        return call_name, result, float(cumulative_t)

    results = PerformanceData()
    total_time = 0
    worst_time = -1

    for gprof_line in output.splitlines():
        gprof_key, gprof_block, cumulative_time = parse_line(gprof_line)

        if gprof_block is None:
            continue

        results.call_list[gprof_key] = gprof_block
        total_time = max(total_time, cumulative_time)

        if gprof_block.perf_time > worst_time:
            worst_time = gprof_block.perf_time
            results.worst_entry = gprof_key, gprof_block

    results.sample_time = total_time * 1000
    results.execution_time = total_time * 1000

    return results


def run_profiler(model: Model, solver: CSolver, trajectories=4, timesteps=101, end_time=100):
    """
    Run low-level performance tests on the specified solver.
    Profiling is performed with a system-level profiling tool.
    Only runtime performance is tested.

    Output is automatically parsed and formatted.
    """
    # The simulation is built without running the solver directly.
    # "Steal" the solver's build engine if it exists, otherwise make one.
    build: BuildEngine = solver.build_engine \
        if isinstance(solver, CSolver) and solver.build_engine is not None \
        else BuildEngine()

    try:
        # build.makefile = MAKEFILE_PATH
        build.prepare(model)

        # Prepare the location where performance data (gmon.out) is written.
        gmon_env = {
            "GMON_OUT_PREFIX": os.path.join(str(build.output_dir), "profile")
        }
        gmon_env.update(os.environ)
        exe = build.build_simulation(simulation_name=solver.target, level=BuildLevel.PROFILE, env=gmon_env)

        # Execute the simulation before running the profiler.
        # When profiling compiler flags (-pg) are enabled,
        #  running the program outputs a profiling data file:
        #  f"{output_dir}/profile.{pid}"
        # This profiler data can then be parsed using gprof.
        process_args = [
            exe,
            "--trajectories", str(trajectories),
            "--timesteps", str(timesteps),
            "--end", str(end_time),
        ]
        start = time.perf_counter()
        subprocess.check_call(" ".join(process_args), env=gmon_env)
        stop = time.perf_counter()

        # Locate gprof profiling metadata in output directory.
        # Gprof will always append `.{pid}` to the end of {GMON_OUT_PREFIX}
        # As such, data is located by finding the file in the output dir
        #   which contains the matching prefix defined above.
        gprof_data: Path = None

        for n in build.output_dir.iterdir():
            if n.stem == "profile":
                gprof_data = build.output_dir.joinpath(n.name)
                break

        if gprof_data is None:
            gprof_data = Path(os.path.join(os.path.dirname(__file__), "gmon.out"))
            if not gprof_data.exists():
                raise EnvironmentError(
                    "Profiler data was not found in the current environment; aborting"
                )

        # Run gprof to process the profiler output.
        # -b = brief, shows minimal output
        # -p = flat profile, discards call graph data
        gprof_args = ["gprof", "-bp", exe, str(gprof_data)]
        perf_out = subprocess.check_output(gprof_args).decode("utf-8")
        gprof_data.unlink()

        # Parse the resulting output
        performance_results = parse_gprof_output(perf_out)
        performance_results.execution_time = (stop - start) * 1000

        return performance_results

    finally:
        build.clean()
