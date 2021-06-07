from gillespy2.core import gillespyError, GillesPySolver, log
from gillespy2.solvers.utilities import solverutils as cutils
import signal  # for solver timeout implementation
import time #for solver timeout implementation
import os #for getting directories for C++ files
import shutil #for deleting/copying files
import threading # for handling read/write to simulation in the background
import subprocess #For calling make and executing c solver
import tempfile #for temporary directories


GILLESPY_PATH = os.path.dirname(os.path.abspath(__file__))
GILLESPY_CPP_TAU_DIR = os.path.join(GILLESPY_PATH, 'c_base/tau_leaping_cpp_solver')
MAKE_FILE = os.path.dirname(os.path.abspath(__file__)) + '/c_base/tau_leaping_cpp_solver/makefile'
CBASE_DIR = os.path.join(GILLESPY_PATH, 'c_base/')


class TauLeapingCSolver(GillesPySolver):
    name = "TauLeapingCSolver"

    def __init__(self, model=None, output_directory=None, delete_directory=True, resume=None, variable=True):
        super(TauLeapingCSolver, self).__init__()
        self.__compiled = False
        self.delete_directory = False
        self.model = model
        self.resume = resume
        self.variable = variable
        if self.model is not None:
            # Create constant, ordered lists for reactions/species/
            self.species_mappings = self.model.sanitized_species_names()
            self.species = list(self.species_mappings.keys())
            self.parameter_mappings = self.model.sanitized_parameter_names()
            self.parameters = list(self.parameter_mappings.keys())
            self.reactions = list(self.model.listOfReactions.keys())

            if isinstance(output_directory, str):
                output_directory = os.path.abspath(output_directory)

                if isinstance(output_directory, str):
                    if not os.path.isfile(output_directory):
                        self.output_directory = output_directory
                        self.delete_directory = delete_directory
                        if not os.path.isdir(output_directory):
                            os.makedirs(self.output_directory)
                    else:
                        raise gillespyError.DirectoryError("File exists with the same path as directory.")
            else:
                self.temporary_directory = tempfile.TemporaryDirectory()
                self.output_directory = self.temporary_directory.name

            if not os.path.isdir(self.output_directory):
                raise gillespyError.DirectoryError("Errors encountered while setting up directory for Solver C++ files."
                                                   )
            self.__write_template()
            self.__compile()

    def __del__(self):
        if self.delete_directory and os.path.isdir(self.output_directory):
            shutil.rmtree(self.output_directory)

    def __write_template(self):
        # Open up template file for reading.

        if self.variable:
            template_file = 'VariableTauSimulationTemplate.cpp'
        else:
            template_file = 'TauSimulationTemplate.cpp'

        # Open up template file for reading.
        with open(os.path.join(GILLESPY_CPP_TAU_DIR, template_file), 'r') as template:
            # Write simulation C++ file.
            template_keyword = "__DEFINE_"
            # Use same lists of model's species and reactions to maintain order
            with open(os.path.join(self.output_directory, 'TauSimulation.cpp'), 'w') as outfile:
                for line in template:
                    if line.startswith(template_keyword):
                        line = line[len(template_keyword):]
                        if line.startswith("VARIABLES"):
                            cutils.write_variables(outfile, self.model, self.reactions, self.species,
                                                   self.parameter_mappings, self.resume, variable=self.variable)
                        if line.startswith("PROPENSITY"):
                            cutils.write_propensity(outfile, self.model, self.species_mappings, self.parameter_mappings
                                                    , self.reactions)
                        if line.startswith("REACTIONS"):
                            cutils.write_reactions(outfile, self.model, self.reactions, self.species)
                        if self.variable:
                            if line.startswith("PARAMETER_UPDATES"):
                                cutils.update_parameters(outfile, self.parameters, self.parameter_mappings)
                    else:
                        outfile.write(line)

    def __compile(self):
        cmd = ["make", "-C", self.output_directory, '-f', MAKE_FILE,
               'TauSimulation', 'GILLESPY_CPP_TAU_DIR=' + GILLESPY_CPP_TAU_DIR, 'CBASE_DIR=' + CBASE_DIR]
        if self.resume:
            if self.resume[0].model != self.model:
                raise gillespyError.ModelError('When resuming, one must not alter the model being resumed.')
        try:
            built = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except KeyboardInterrupt:
            log.warning("Solver has been interrupted during compile time, unexpected behavior may occur.")
        if built.returncode == 0:
            self.__compiled = True
        else:
            raise gillespyError.BuildError("Error encountered while compiling file:\nReturn code: "
                                           "{0}.\nError:\n{1}\n{2}\n".format(built.returncode,
                                                                             built.stdout.decode('utf-8'),
                                                                             built.stderr.decode('utf-8')))

    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'timeout', 'increment', 'seed', 'debug', 'profile')

    def run(self=None, model=None, t=20, number_of_trajectories=1, timeout=0,
            increment=0.05, seed=None, debug=False, profile=False, resume=None,
            tau_step=.03, variables={}, tau_tol=0.03, **kwargs):

        pause = False
        if resume is not None:
            if t < resume['time'][-1]:
                raise gillespyError.ExecutionError(
                    "'t' must be greater than previous simulations end time, or set in the run() method as the "
                    "simulations next end time")

        if self is None or self.model is None:
            self = TauLeapingCSolver(model, resume=resume)


        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        unsupported_sbml_features = {
            'Rate Rules': len(model.listOfRateRules),
            'Assignment Rules': len(model.listOfAssignmentRules),
            'Events': len(model.listOfEvents),
            'Function Definitions': len(model.listOfFunctionDefinitions)
        }
        detected_features = []
        for feature, count in unsupported_sbml_features.items():
            if count:
                detected_features.append(feature)

        if len(detected_features):
                raise gillespyError.ModelError(
                'Could not run Model.  SBML Feature: {} not supported by TauLeapingSolver.'.format(detected_features))

        if not isinstance(variables, dict):
            raise gillespyError.SimulationError(
                'argument to variables must be a dictionary.')
        for v in variables.keys():
            if v not in self.species+self.parameters:
                raise gillespyError.SimulationError('Argument to variable "{}" \
                is not a valid variable.  Variables must be model species or parameters.'.format(v))

        if self.__compiled:
            self.simulation_data = None
            if resume is not None:
                t = abs(t - resume['time'][-1])

            number_timesteps = int(round(t / increment + 1))

            # Execute simulation.
            args = [os.path.join(self.output_directory, 'TauSimulation'), '-trajectories', str(number_of_trajectories),
                    '-timesteps', str(number_timesteps), '-tau_step', str(tau_step), '-end', str(t),
                    '-tau_tol', str(tau_tol)]

            if self.variable:  # Is a variable simulation
                populations = cutils.update_species_init_values(model.listOfSpecies, self.species, variables, resume)
                parameter_values = cutils.change_param_values(model.listOfParameters, self.parameters, model.volume,
                                                              variables)
                args.extend(['-initial_values', populations, '-parameters', parameter_values])

            if seed is not None:
                if isinstance(seed, int):
                    args.append('-seed')
                    args.append(str(seed))
                else:
                    seed_int = int(seed)
                    if seed_int > 0:
                        args.append('-seed')
                        args.append(str(seed_int))
                    else:
                        raise gillespyError.ModelError("seed must be a positive integer")

            # Handler for reading data from subprocess, in background thread.
            def sim_delegate(sim, sim_buffer):
                def read_next():
                    # Reads the next block from the simulation output.
                    # Returns the length of the string read.
                    line = sim.stdout.read().decode("utf-8")
                    ln = len(line)
                    if ln > 0:
                        sim_buffer.append(line)
                    return ln

                # Read output 1 block at a time, until the program is finished.
                page_size = read_next()
                while page_size > 0 and sim.poll() is None:
                    page_size = read_next()

            # Buffer to store the output of the simulation (retrieved from sim_delegate thread).
            buffer = []
            # Each platform is given their own platform-specific sub_kill() function.
            # Windows event handling
            if os.name == "nt":
                sub_kill = lambda sim: sim.send_signal(signal.CTRL_BREAK_EVENT)
                platform_args = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP,
                                 "start_new_session": True}
            # POSIX event handling
            else:
                sub_kill = lambda sim: os.killpg(sim.pid, signal.SIGINT)
                platform_args = {"start_new_session": True}

            thread_events = {"timeout": False}
            with subprocess.Popen(args, stdout=subprocess.PIPE, **platform_args) as simulation:
                # Put a timer on in the background, if a timeout was specified.
                def timeout_kill():
                    thread_events["timeout"] = True
                    sub_kill(simulation)

                timeout_thread = threading.Timer(timeout, timeout_kill)
                try:
                    output_process = threading.Thread(name="SimulationHandlerThread",
                                                      target=sim_delegate,
                                                      args=(simulation, buffer))
                    output_process.start()
                    if timeout > 0:
                        timeout_thread.start()

                    # Poll for the program's status; keyboard interrupt is ignored if we use .wait()
                    while simulation.poll() is None:
                        time.sleep(0.1)
                except KeyboardInterrupt:
                    sub_kill(simulation)
                    pause = True
                finally:
                    # Finish off the output reader thread and the timer thread.
                    return_code = simulation.wait()
                    output_process.join()
                    if timeout_thread.is_alive():
                        timeout_thread.cancel()

                    # Decode from byte, split by comma into array
                    stdout = "".join(buffer).split(",")
                    # Check if the simulation had been paused
                    # (Necessary because we can't set pause to True from thread handler)
                    if thread_events["timeout"]:
                        pause = True
                        return_code = 33

            if return_code in [0, 33]:
                trajectory_base, timeStopped = cutils.parse_binary_output(number_of_trajectories,
                                                                          number_timesteps,
                                                                          len(model.listOfSpecies), stdout,
                                                                          pause=pause)
                if model.tspan[2] - model.tspan[1] == 1:
                    timeStopped = int(timeStopped)

                # Format results
                self.simulation_data = []
                for trajectory in range(number_of_trajectories):
                    data = {'time': trajectory_base[trajectory, :, 0]}
                    for i in range(len(self.species)):
                        data[self.species[i]] = trajectory_base[trajectory, :, i + 1]

                    self.simulation_data.append(data)
            else:
                raise gillespyError.ExecutionError("Error encountered while running simulation C++ file:"
                                                   "\nReturn code: {0}.\nError:\n{1}\n".
                                                   format(simulation.returncode, simulation.stderr))

            if resume is not None or timeStopped != 0:
                self.simulation_data = cutils.c_solver_resume(timeStopped, self.simulation_data, t,
                                                              resume=resume)

        return self.simulation_data, return_code