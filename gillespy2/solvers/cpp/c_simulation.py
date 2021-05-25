import os
import subprocess
import signal
import threading

import numpy

from enum import IntEnum
from concurrent.futures import Future, ThreadPoolExecutor

from gillespy2.core import Model
from gillespy2.core import gillespyError
from gillespy2.solvers.cpp.c_decoder import SimDecoder
from gillespy2.solvers.cpp.build.build_engine import BuildEngine

class SimulationReturnCode(IntEnum):
    DONE = 0
    PAUSED = 33
    FAILED = -1

class CSimulation:
    def __init__(self, model=None, output_directory=None, delete_directory=True, resume=None, variable=True):
        self.delete_directory = False
        self.model = model
        self.resume = resume
        self.variable = variable

        # Validate output_directory, ensure that it doesn't already exist
        if isinstance(output_directory, str):
            output_directory = os.path.abspath(output_directory)

            if os.path.exists(output_directory):
                raise gillespyError.DirectoryError(
                    f"Could not write to specified output directory: {output_directory}"
                    " (already exists)"
                )
            self.output_directory = output_directory
            self.delete_directory = delete_directory

        if self.model is None:
            return

        self.species_mappings = self.model.sanitized_species_names()
        self.species = list(self.species_mappings.keys())
        self.parameter_mappings = self.model.sanitized_parameter_names()
        self.parameters = list(self.parameter_mappings.keys())
        self.reactions = list(self.model.listOfReactions.keys())

    def _build(self, model: Model, simulation_name: str, variable: bool, debug: bool = False) -> str:
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

    def _run_async(self, sim_exec: str, sim_args: list[str], decoder: SimDecoder, timeout: int = 0) -> Future[int]:
        """
        Run the target executable simulation as async.

        :param sim_exec: The executable simulation to run.
        :type sim_exec: str

        :returns: A future which represents the currently executing run_simulation job.
        """

        executor = ThreadPoolExecutor()
        return executor.submit(self._run, sim_exec, sim_args, decoder, timeout)

    def _run(self, sim_exec: str, sim_args: list[str], decoder: SimDecoder, timeout: int = 0) -> int:
        """
        Run the target executable simulation.

        :param sim_exec: The executable simulation to run.
        :type sim_exec: str

        :returns: The return_code of the simulation.
        """

        # Prefix the executable to the sim arguments.
        sim_args = [sim_exec] + sim_args
        print(sim_args)

        # nt and *nix require different methods to force shutdown a running process.
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

        timeout_event = [False]

        with subprocess.Popen(sim_args, stdout=subprocess.PIPE, **platform_args) as simulation:
            def timeout_kill():
                timeout_event[0] = True
                proc_kill(simulation)

            timeout_thread = threading.Timer(timeout, timeout_kill)

            if timeout > 0:
                timeout_thread.start()

            try:
                decoder.read(simulation.stdout)
            
            except KeyboardInterrupt:
                proc_kill(simulation)

            finally:
                return_code = simulation.wait()

                if timeout_thread.is_alive():
                    timeout_thread.cancel()

                if timeout_event[0]:
                    return SimulationReturnCode.PAUSED

                print(return_code)

                if return_code not in [0, 33]:
                    return SimulationReturnCode.FAILED

                return SimulationReturnCode.DONE

    def _make_args(self, args_dict: dict[str, str]) -> list[str]:
        args_list = []

        for key, value in args_dict.items():
            args_list.extend([f"-{key}", str(value)])

        return args_list

    def _format_output(self, trajectories: numpy.ndarray):
        # Check the dimensionality of the input trajectory.
        if not len(trajectories.shape) == 3:
            raise gillespyError.ValidationError("Could not format trajectories, input numpy.ndarray is not 3-dimensional.")

        # The trajectory count is the first dimention of the input ndarray.
        trajectory_count = trajectories.shape[0]
        simulation_data = [ ]

        # Begin iterating through the trajectories, copying each dimension into simulation_data.
        for trajectory in range(trajectory_count):
            # Copy the first index of the third-dimension into the simulation_data dictionary.
            data = {
                "time": trajectories[trajectory, :, 0]
            }

            for i in range(len(self.species)):
                data[self.species[i]] = trajectories[trajectory, :, i + 1]

            simulation_data.append(data)

        return simulation_data
