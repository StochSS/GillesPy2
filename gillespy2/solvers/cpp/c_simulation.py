from enum import Enum
import os
import subprocess
import signal
import threading

from enum import Enum
from concurrent.futures import Future, ThreadPoolExecutor
from turtle import done
from gillespy2.solvers.cpp.c_encoder import SimDecoder

from gillespy2.core import Model
from gillespy2.solvers.cpp.build.build_engine import BuildEngine

class SimulationReturnCode(Enum):
    DONE = 1
    PAUSED = 2
    FAILED = 3

def build(model: Model, simulation_name: str, variable: bool, debug: bool = False) -> str:
    """
    Generate and build the simulation from the specified Model and solver_name into the output_dir.

    :param model: The Model to simulate.
    :type model: gillespy2.Model

    :param solver_name: The name of the simulation to execute.
    :type str:

    :param output_directory: The directory to output the simulation executable.
    :type str:
    """

    # Prepare the build workspace.
    build_engine = BuildEngine(debug=debug)
    build_engine.prepare(model, variable)

    # Compile the simulation, returning the path of the executable.
    return build_engine.build_simulation(simulation_name)

def run_async(sim_exec: str, sim_args: list[str], decoder: SimDecoder, timeout: int = 0) -> Future[int]:
    """
    Run the target executable simulation as async.

    :param sim_exec: The executable simulation to run.
    :type sim_exec: str

    :returns: A future which represents the currently executing run_simulation job.
    """

    executor = ThreadPoolExecutor()
    return executor.submit(run, sim_exec, sim_args, decoder)

def run(sim_exec: str, sim_args: list[str], decoder: SimDecoder, timeout: int = 0) -> int:
    """
    Run the target executable simulation.

    :param sim_exec: The executable simulation to run.
    :type sim_exec: str

    :returns: The return_code of the simulation.
    """

    # Prefix the executable to the sim arguments.
    sim_args = [sim_exec] + sim_args

    if os.name == "nt":
        proc_kill = lambda sim: sim.send_signal(signal.CTRL_BREAK_EVENT)
        platform_args = {
            "creationargs": subprocess.CREATE_NEW_PROCESS_GROUP,
            "start_new_session": True
        }

    else:
        proc_kill = lambda sim: os.killpg(sim.pid, signal.SIGINT)
        platform_args = {
            "start_new_session": True
        }

    pause = False
    timeout_event = [False]

    with subprocess.Popen(sim_args, stdout=subprocess.PIPE, **platform_args) as simulation:
        def timeout_kill():
            timeout_event[0] = True
            proc_kill(simulation)

        if timeout > 0:
            timeout_thread = threading.Timer(timeout, timeout_kill)
            timeout_thread.start()

        try:
            decoder.read(simulation.stdout)
        
        except KeyboardInterrupt:
            proc_kill(simulation)
            pause = True

        finally:
            return_code = simulation.wait()

            if timeout_thread.is_alive():
                timeout_thread.cancel()

            if timeout_event[0]:
                return SimulationReturnCode.PAUSED

            if return_code not in [0, 33]:
                return SimulationReturnCode.FAILED

            return SimulationReturnCode.DONE
