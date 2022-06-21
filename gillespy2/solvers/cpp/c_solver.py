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
import copy
import subprocess
import signal
import threading
import queue
import platform

import numpy

from enum import IntEnum
from concurrent.futures import ThreadPoolExecutor
from typing import Union

import gillespy2
from gillespy2.core import Model
from gillespy2.core import Results
from gillespy2.core import log
from gillespy2.core import gillespyError
from gillespy2.solvers.cpp.c_decoder import SimDecoder
from gillespy2.solvers.cpp.build.build_engine import BuildEngine
from gillespy2.solvers.cpp.build.template_gen import SanitizedModel

class SimulationReturnCode(IntEnum):
    DONE = 0
    PAUSED = 33
    FAILED = -1

class CSolver:
    """
    This class implements base behavior that will be needed for C++ solver implementees.

    :param model: The Model to simulate.
    :type model: Model

    :param output_directory: The working output directory.
    :type output_directory: str

    :param delete_directory: If True then the output_directory will be deleted upon completetion.
    :type delete_directory: bool

    :param resume: Resume data from a previous simulation run.

    :param variable: Indicates whether the simulation should be variable.
    :type variable: bool
    """
    rc = 0

    def __init__(self, model: Model = None, output_directory: str = None, delete_directory: bool = True, resume=None, variable: bool = False):
        if model is None:
            raise gillespyError.SimulationError("A model is required to run the simulation.")

        if len(BuildEngine.get_missing_dependencies()) > 0:
            raise gillespyError.SimulationError(
                "Please install/configure 'g++' and 'make' on your system, to ensure that GillesPy2 C solvers will run properly."
            )

        if platform.system() == "Windows" and " " in gillespy2.__file__:
            raise gillespy2Error.SimulationError("GillesPy2 does not support spaces in its path on Windows systems.")

        self.delete_directory = False
        self.model = copy.deepcopy(model)
        self.resume = resume
        self.variable = variable
        self.build_engine: BuildEngine = None

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

        if self.model is not None:
            self._set_model()

        self.is_instantiated = True

    def __del__(self):
        if self.build_engine is None:
            return

        if not self.delete_directory:
            return

        self.build_engine.clean()

    def _build(self, model: "Union[Model, SanitizedModel]", simulation_name: str, variable: bool, debug: bool = False) -> str:
        """
        Generate and build the simulation from the specified Model and solver_name into the output_dir.

        :param model: The Model to simulate.
        :type model: gillespy2.Model

        :param simulation_name: The name of the simulation to execute.
        :type simulation_name: str

        :param variable: If True the simulation will be variable, False if not.
        :type variable: bool

        :param debug: Enables or disables debug behavior.
        :type debug: bool
        """

        # Prepare the build workspace.
        if self.build_engine is None or self.build_engine.get_executable_path() is None:
            self.build_engine = BuildEngine(debug=debug, output_dir=self.output_directory)
            self.build_engine.prepare(model, variable)
            # Compile the simulation, returning the path of the executable.
            return self.build_engine.build_simulation(simulation_name)

        # Assume that the simulation has already been built.
        return self.build_engine.get_executable_path()

    def _run_async(self, sim_exec: str, sim_args: "list[str]", decoder: SimDecoder, timeout: int = 0):
        """
        Run the target executable simulation as async.

        :param sim_exec: The executable simulation to run.
        :type sim_exec: str

        :param sim_args: The arguments to pass on simulation run.
        :type sim_args: list[str]

        :param decoder: The SimDecoder instance that will handle simulation output.
        :type decoder: SimDecoder

        :returns: A future which represents the currently executing run_simulation job.
        """

        executor = ThreadPoolExecutor()
        return executor.submit(self._run, sim_exec, sim_args, decoder, timeout)

    def _run(self, sim_exec: str, sim_args: "list[str]", decoder: SimDecoder, timeout: int = 0, display_args: dict = None) -> int:
        """
        Run the target executable simulation.

        :param sim_exec: The executable simulation to run.
        :type sim_exec: str

        :param sim_args: The arguments to pass on simulation run.
        :type sim_args: list[str]

        :param decoder: The SimDecoder instance that will handle simulation output.
        :type decoder: SimDecoder

        :param display_args: The kwargs need to setup the live graphing
        :type display_args: dict

        :returns: The return_code of the simulation.
        """

        # Prefix the executable to the sim arguments.
        sim_args = [sim_exec] + sim_args

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

        live_grapher = [None]

        if display_args is not None:
            live_queue = queue.Queue(maxsize=1)
            def decoder_cb(curr_time, curr_state, trajectory_base=decoder.trajectories):
                try:
                    old_entry = live_queue.get_nowait()
                except queue.Empty as err:
                    pass
                curr_state = {self.species[i]: curr_state[i] for i in range(len(curr_state))}
                entry = ([curr_state], [curr_time], trajectory_base)
                live_queue.put(entry)
                
            decoder.with_callback(decoder_cb)

            from gillespy2.core.liveGraphing import (
                LiveDisplayer, CRepeatTimer, valid_graph_params
            )
            valid_graph_params(display_args['live_output_options'])
            live_grapher[0] = LiveDisplayer(**display_args)
            display_timer = CRepeatTimer(
                display_args['live_output_options']['interval'], live_grapher[0].display,
                args=(live_queue, display_args['live_output_options']['type'])
            )

        timeout_event = [False]
        with subprocess.Popen(sim_args, stdout=subprocess.PIPE, **platform_args) as simulation:
            def timeout_kill():
                timeout_event[0] = True
                proc_kill(simulation)

            timeout_thread = threading.Timer(timeout, timeout_kill)
            reader_thread = threading.Thread(target=decoder.read,
                                             args=(simulation.stdout,))

            if timeout > 0:
                timeout_thread.start()

            try:
                reader_thread.start()
                if display_args is not None:
                    display_timer.start()
                reader_thread.join()
            
            except KeyboardInterrupt:
                if live_grapher[0] is not None:
                    display_timer.pause = True
                proc_kill(simulation)

            finally:
                if live_grapher[0] is not None:
                    display_timer.cancel()
                timeout_thread.cancel()
                return_code = simulation.wait()
                reader_thread.join()

                if timeout_event[0]:
                    return SimulationReturnCode.PAUSED

                return self._handle_return_code(return_code)

    def _make_args(self, args_dict: "dict[str, str]") -> "list[str]":
        """
        Convert a dictionary of key, value pairs into a valid Popen argument list.
        Note: Do not prefix a key with `-` as this will be handled automatically.

        :param args_dict: A dictionary of named arguments.
        :type args_dict: dict[str, str]

        :returns: A formatted list of arguments.
        """

        args_list = []

        for key, value in args_dict.items():
            args_list.extend([f"--{key}", str(value)])

        return args_list

    def _format_output(self, trajectories: numpy.ndarray):
        # Check the dimensionality of the input trajectory.
        if not len(trajectories.shape) == 3:
            raise gillespyError.ValidationError("Could not format trajectories, input numpy.ndarray is not 3-dimensional.")

        # The trajectory count is the first dimention of the input ndarray.
        trajectory_count = trajectories.shape[0]
        self.result = []

        # Begin iterating through the trajectories, copying each dimension into simulation_data.
        for trajectory in range(trajectory_count):
            # Copy the first index of the third-dimension into the simulation_data dictionary.
            data = {
                "time": trajectories[trajectory, :, 0]
            }

            for i in range(len(self.species)):
                data[self.species[i]] = trajectories[trajectory, :, i + 1]

            self.result.append(data)

        return self.result

    def _handle_return_code(self, return_code: "int") -> "SimulationReturnCode":
        """
        Default return code handler; determines whether the simulation succeeded or failed.
        Intended to be overridden by solver subclasses, which handles solver-specific return codes.

        Does nothing if the return code checks out, otherwise raises an error.

        :param return_code: Return code returned by a simulation.
        :type return_code: int
        """
        if return_code == 33:
            return SimulationReturnCode.PAUSED
        if return_code == 0:
            return SimulationReturnCode.DONE

        raise gillespyError.ExecutionError("Error encountered while running simulation C++ file "
                                           f"(return code: {int(return_code)})")

    def _make_resume_data(self, time_stopped: int, simulation_data: numpy.ndarray, t: int):
        """
        If the simulation was paused then the output data needs to be trimmed to allow for resume.
        In the event the simulation was not paused, no data is changed.
        """
        # No need to create resume data if the simulation completed without interruption.
        # Note that, currently, some C++ solvers write 0 out as the "time stopped" by default.
        # This is likely to change in the future.
        if not time_stopped < t or time_stopped == 0:
            return simulation_data

        # Find the index of the time step value which is closest to the time stopped.
        cutoff = numpy.searchsorted(simulation_data[-1]["time"], float(time_stopped))
        if cutoff < 2:
            log.warning('You have paused the simulation too early, and no points have been calculated past'
                        ' initial values. A graphic display will not produce expected results.')

        # Break off any extraneous data which goes past the cutoff time.
        # Any data in this case is assumed to be untrusted.
        for entry_name, entry_data in simulation_data[-1].items():
            simulation_data[-1][entry_name] = entry_data[:cutoff]

        return simulation_data

    def _set_model(self, model=None):
        if model is not None:
            self.model = copy.deepcopy(model)

        self._build(self.model, self.target, self.variable, False)
        self.species_mappings = self.model.sanitized_species_names()
        self.species = list(self.species_mappings.keys())
        self.parameter_mappings = self.model.sanitized_parameter_names()
        self.parameters = list(self.parameter_mappings.keys())
        self.reactions = list(self.model.listOfReactions.keys())
        self.result = []
        self.rc = 0

    def _update_resume_data(self, resume: Results, simulation_data: "list[dict[str, numpy.ndarray]]", time_stopped: int):
        """
        Modify the simulation output to continue from a previous Results object.
        Does not handle the case where the simulation was interrupted again.
        """
        # No need to update the resume data if there is no previous data, or not enough.
        if resume is None or len(resume["time"]) < 2:
            return simulation_data

        resume_time = float(resume["time"][-1])
        increment = resume_time - float(resume["time"][-2])
        # Replace the simulation's timespan to continue where the Results object left off.
        simulation_data[-1]["time"] = numpy.arange(start=(resume_time),
                                                   stop=(resume_time + time_stopped + increment),
                                                   step=increment)

        for entry_name, entry_data in simulation_data[-1].items():
            # The results of the current simulation is treated as an "extension" of the resume data.
            # As such, the new simulation output is formed by joining the two end to end.
            new_data = numpy.concatenate((resume[entry_name], entry_data[1:]), axis=None)
            simulation_data[-1][entry_name] = new_data

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
            log.warning(f"Unsupported keyword argument for solver {self.name}: {key}")

    def _validate_seed(self, seed: int):
        if seed is None:
            return None

        if not isinstance(seed, int):
            seed = int(seed)

        if seed <= 0:
            raise gillespyError.ModelError("`seed` must be a postive integer.")

        return seed
        
    def _validate_variables_in_set(self, variables, values):
        for var in variables.keys():
            if var not in values:
                raise gillespyError.SimulationError(f"Argument to variable '{var}' is not a valid variable. Variables must be model species or parameters.")

    def _validate_type(self, value, typeof: type, message: str):
        if not type(value) == typeof:
            raise gillespyError.SimulationError(message)
