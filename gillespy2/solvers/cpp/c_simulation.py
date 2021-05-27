from gillespy2.solvers.utilities import solverutils
import logging
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
    def __init__(self, model: Model = None, output_directory: str = None, delete_directory: bool = True, resume=None, variable: bool = False):
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
        self.build_engine = None

        if self.model is None:
            return

        self.species_mappings = self.model.sanitized_species_names()
        self.species = list(self.species_mappings.keys())
        self.parameter_mappings = self.model.sanitized_parameter_names()
        self.parameters = list(self.parameter_mappings.keys())
        self.reactions = list(self.model.listOfReactions.keys())
        self.simulation_data = []

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
        self.build_engine = BuildEngine(debug=debug, output_dir=self.output_directory)
        self.build_engine.prepare(model, variable)

        # Compile the simulation, returning the path of the executable.
        return self.build_engine.build_simulation(simulation_name)

    def _run_async(self, sim_exec: str, sim_args: "list[str]", decoder: SimDecoder, timeout: int = 0) -> "Future[int]":
        """
        Run the target executable simulation as async.

        :param sim_exec: The executable simulation to run.
        :type sim_exec: str

        :returns: A future which represents the currently executing run_simulation job.
        """

        executor = ThreadPoolExecutor()
        return executor.submit(self._run, sim_exec, sim_args, decoder, timeout)

    def _run(self, sim_exec: str, sim_args: "list[str]", decoder: SimDecoder, timeout: int = 0) -> int:
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
                "creationflags": subprocess.CREATE_NEW_PROCESS_GROUP,
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

                # Clean up if specified to do so by the user
                if self.build_engine is not None and self.delete_directory:
                    self.build_engine.clean()

                if return_code not in [0, 33]:
                    return SimulationReturnCode.FAILED

                return SimulationReturnCode.DONE

    def _make_args(self, args_dict: "dict[str, str]") -> "list[str]":
        """
        Convert a dictionary of key, value pairs into a valid Popen argument list.
        Note: Do not prefix a key with `-` as this will be handled automatically.

        :param args_dict: A dictionary of named arguments.

        :returns: A formatted list of arguments.
        """

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
        self.simulation_data = []

        # Begin iterating through the trajectories, copying each dimension into simulation_data.
        for trajectory in range(trajectory_count):
            # Copy the first index of the third-dimension into the simulation_data dictionary.
            data = {
                "time": trajectories[trajectory, :, 0]
            }

            for i in range(len(self.species)):
                data[self.species[i]] = trajectories[trajectory, :, i + 1]

            self.simulation_data.append(data)

        return self.simulation_data

    def _make_resume_data(self, time_stopped: int, simulation_data: numpy.ndarray, t: int, resume):
        """
        If the simulation was paused then the output data needs to be trimmed to allow for resume.
        In the event the simulation was not paused, no data is changed.
        """

        if resume is not None or time_stopped != 0:
            return solverutils.c_solver_resume(time_stopped, simulation_data, t, resume=resume)

        return simulation_data

    def _validate_resume(self, t: int, resume):
        """
        Validate `resume`. An exception will be thrown if resume['time'][-1] is > t.
        """

        if resume is None:
            return

        if t < resume["time"][-1]:
            raise gillespyError.ExecutionError(
                "'t' must be greater than previous simulations end time, or set in the run() method as the "
                "simulations next end time"
            )

    def _validate_kwargs(self, **kwargs):
        """
        Validate any additional kwargs passed to the model. If any exist, warn the user.
        """

        if len(kwargs) == 0:
            return

        for key, val in kwargs.items():
            logging.warn(f"Unsupported keyword argument for solver {self.name}: {key}")

    def _validate_sbml_features(self, unsupported_features: "dict[str, str]"):
        detected = [ ]
        for feature_name, count in unsupported_features.items():
            if count:
                detected.append(feature_name)

        if len(detected):
            raise gillespyError.ModelError(f"Could not run Model.  SBML Feature: {detected} not supported by SSACSolver.")

    def _validate_seed(self, seed: int):
        if seed is None:
            return None

        if not isinstance(seed, int):
            seed = int(seed)

        if seed <= 0:
            raise gillespyError.ModelError("`seed` must be a postive integer.")

        return seed
        
    def _validate_variables_in_set(self, variables, set):
        for var in variables.keys():
            if var not in set:
                raise gillespyError.SimulationErrorp(f"Argument to variable '{var}' is not a valid variable. Variables must be model species or parameters.")

    def _validate_type(self, value, typeof: type, message: str):
        if not type(value) == typeof:
            raise gillespyError.SimulationError(message)
