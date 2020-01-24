import numpy as np
import random
import tempfile
import os
import uuid
import subprocess
import shutil
from gillespy2.core import GillesPySolver, Model, log
from gillespy2.core.gillespyError import SimulationError, InvalidModelError


class StochKitBaseSolver(GillesPySolver):
    name = "StochKitBaseSolver"
    """
    Abstract class for a StochKit solver. This is generally called from within a
    GillesPy Model through the Model.run function. Returns simulation
    trajectories.

    Attributes
    ----------
    model : gillespy.Model
        The model on which the solver will operate.
    t : float
        The end time of the solver.
    number_of_trajectories : int
        The number of times to sample the chemical master equation. Each
        trajectory will be returned at the end of the simulation.
    increment : float
        The time step of the solution.
    seed : int
        The random seed for the simulation. Defaults to None.
    stochkit_home : str
        Path to stochkit. This is set automatically upon installation, but
        may be overwritten if desired.
    algorithm : str
        The solver by which to simulate the model. 'ssa' or 'tau_leaping'
        are the available options. If 'ssa' is chosen, StochKit will choose
        from the available ssa options.
    job_id : str
        If given, this will be the name of the solver run. Usually not set.
    extra_args : str
        Any extra arguments for the stochkit solver. See StochKit2
        documentation for details.
    debug : bool (False)
        Set to True to provide additional debug information about the
        simulation.
    show_labels : bool (True)
        Use names of species as index of result object rather than position numbers.
    """
    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
            stochkit_home=None, algorithm=None, job_id=None, extra_args='',
            debug=False, profile=False, show_labels=False, **kwargs):
        """
        Call out and run the solver. Collect the results.
        """
        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to solver: {0}'.format(key))

        if algorithm is None:
            raise SimulationError("No algorithm selected")

        # We write all StochKit input and output files to a temporary folder
        prefix_base_dir = tempfile.mkdtemp()
        prefix_out_dir = os.path.join(prefix_base_dir, 'output')
        os.mkdir(prefix_out_dir)

        if job_id is None:
            job_id = str(uuid.uuid4())

        if isinstance(model, Model):
            # Write a temporary StochKit2 input file.
            outfile = os.path.join(prefix_base_dir, "temp_input_{0}.xml".format(job_id))
            with open(outfile, 'w') as model_file_handle:
                model_file_handle.write(model.serialize())
        elif isinstance(model, str):
            outfile = model
        else:
            raise InvalidModelError('Model must be either a GillesPy Model instance or an xml file name.')

        executable = cls.locate_executable(stochkit_home=stochkit_home, algorithm=algorithm)

        if executable is None:
            raise SimulationError("stochkit executable '{0}' not found. \
                Make sure it is your path, or set STOCHKIT_HOME environment \
                variable'".format(algorithm))

        # Assemble argument list for StochKit
        out_dir = os.path.join(prefix_out_dir, job_id)

        if increment is None:
            increment = t / 20.0
        num_output_points = round(t / increment)

        # Assemble the argument list
        args = '--model {0} --out-dir {1} -t {2} -i {3}'.format(outfile, out_dir, t, int(num_output_points))

        directories = os.listdir(prefix_out_dir)
        if os.path.isdir(out_dir):
            if debug:
                print('Ensemble {0} already existed, using --force.'.format(job_id))
            args += ' --force'

        # If we are using local mode, shell out and run StochKit
        # (SSA or Tau-leaping or ODE)
        cmd = ' '.join([executable, args, extra_args])
        if debug:
            print("cmd: {0}".format(cmd))

        # Execute
        try:
            handle = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            return_code = handle.wait()
        except OSError as e:
            raise SimulationError("Solver execution failed: {0}\n{1}".format(cmd, e))

        try:
            stderr = handle.stderr.read()
        except Exception as e:
            stderr = 'Error reading stderr: {0}'.format(e)
        try:
            stdout = handle.stdout.read()
        except Exception as e:
            stdout = 'Error reading stdout: {0}'.format(e)

        if return_code != 0:
            raise SimulationError("Solver execution failed: '{0}' output: {1}{2}".format(cmd, stdout, stderr))

        try:
            # Get data using solver specific function
            trajectories = cls.get_trajectories(out_dir, debug=debug, show_labels=show_labels)
            if len(trajectories) == 0:
                raise SimulationError("Solver execution failed: '{0}' output: {1}{2}".format(cmd, stdout, stderr))
            if show_labels:
                labels, trajectories = trajectories
                trajectories = cls.label_trajectories(trajectories, labels)
            return trajectories
        except Exception as e:
            compile_log_file = os.path.join(prefix_base_dir, 'temp_input_{0}_generated_code'.format(job_id),
                                            'compile-log.txt')
            log_file = os.path.join(prefix_out_dir, job_id, 'log.txt')
            for file_name in [compile_log_file, log_file]:
                if os.path.isfile(file_name):
                    with open(file_name) as f:
                        error = f.read()
                    raise SimulationError("Error running simulation: {0}\n{1}\n".format(file_name, error))

            raise SimulationError("Error using solver.get_trajectories('{0}'): {1}".format(out_dir, e))
        finally:
            # Clean up
            if debug:
                print("prefix_base_dir={0}".format(prefix_base_dir))
                print("STDOUT: {0}".format(stdout))
                print("STDERR: {0}".format(stderr))
            else:
                shutil.rmtree(prefix_base_dir)

    @staticmethod
    def locate_executable(stochkit_home=None, algorithm=None):
        # Algorithm, SSA or Tau-leaping?
        executable = None
        if stochkit_home is not None:
            if os.path.isfile(os.path.join(stochkit_home, algorithm)):
                executable = os.path.join(stochkit_home, algorithm)
            else:
                raise SimulationError("stochkit executable '{0}' not found \
                stochkit_home={1}".format(algorithm, stochkit_home))
        elif os.environ.get('STOCHKIT_HOME') is not None:
            if os.path.isfile(os.path.join(os.environ.get('STOCHKIT_HOME'),
                                           algorithm)):
                executable = os.path.join(os.environ.get('STOCHKIT_HOME'),
                                          algorithm)
        if executable is None:
            # try to find the executable in the path
            if os.environ.get('PATH') is not None:
                for directory in os.environ.get('PATH').split(os.pathsep):
                    if os.path.isfile(os.path.join(directory, algorithm)):
                        executable = os.path.join(directory, algorithm)
                        break
        return executable

    @staticmethod
    def process_seed(seed):
        if seed is None:
            seed = random.randint(0, 2147483647)
        # StochKit breaks for long ints
        if seed.bit_length() >= 32:
            seed = seed & ((1 << 32) - 1)
            if seed > (1 << 31) - 1:
                seed -= 1 << 32
        return seed

    @staticmethod
    def label_trajectories(trajectories, labels):
        results = []
        for r in trajectories:
            ret = {}
            for n, l in enumerate(labels):
                ret[l] = r[:, n]
            results.append(ret)
        return results


class StochKitSolver(StochKitBaseSolver):
    name = 'StochKitSolver'
    """
    Abstract class for StochKit solver derived from the GillesPySolver class.
    This is generally used to set up the solver.

    Attributes
    ----------
    model : gillespy.Model
        The model on which the solver will operate.
    t : float
        The end time of the solver.
    number_of_trajectories : int
        The number of times to sample the chemical master equation. Each
        trajectory will be returned at the end of the simulation.
    increment : float
        The time step of the solution.
    seed : int
        The random seed for the simulation. Defaults to None.
    stochkit_home : str
        Path to stochkit. This is set automatically upon installation, but
        may be overwritten if desired.
    algorithm : str
        The solver by which to simulate the model. 'ssa' or 'tau_leaping'
        are the available options. If 'ssa' is chosen, StochKit will choose
        from the available ssa options.
    job_id : str
        If given, this will be the name of the solver run. Usually not set.
    method : str
        The specific SSA to call. NOT YET FUNCTIONAL.
    debug : bool (False)
        Set to True to provide additional debug information about the
        simulation.
    """
    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
            stochkit_home=None, algorithm='ssa', job_id=None, method=None,
            debug=False, show_labels=False, profile=False, processes=1, **kwargs):

        # all this is specific to StochKit
        if model.units == "concentration":
            raise SimulationError("StochKit can only simulate population \
                                  models, please convert to population-based model for \
                                  stochastic simulation. Use solver = StochKitODESolver \
                                  instead to simulate a concentration model deterministically.")

        seed = super().process_seed(seed)

        # We keep all the trajectories by default.
        args = ' -p {0} --keep-trajectories --label --seed {1} --realizations {2}'.format(processes, seed, number_of_trajectories)

        if method is not None:  # This only works for StochKit 2.1
            args += ' --method ' + str(method)

        return super().run(model=model, t=t, number_of_trajectories=number_of_trajectories, increment=increment, seed=seed, stochkit_home=stochkit_home,
                           algorithm=algorithm, job_id=job_id, debug=debug, show_labels=show_labels, extra_args=args)

    @classmethod
    def get_trajectories(cls, out_dir, debug=False, show_labels=False):
        if debug:
            print("StochKitSolver.get_trajectories(out_dir={0}".format(out_dir))
        # Collect all the output data
        trajectories = []
        trajectory_directory = os.path.join(out_dir, 'trajectories')
        for filename in os.listdir(trajectory_directory):
            if 'trajectory' in filename:
                filename = os.path.join(trajectory_directory, filename)
                trajectories.append(np.loadtxt(filename, skiprows=1))
            else:
                raise SimulationError("Couldn't identify file '{0}' found in \
                                        output folder".format(filename))
        if show_labels:
            with open(os.path.join(trajectory_directory, 'trajectory0.txt'), 'r') as fd:
                headers = fd.readline()
            return headers.split(), trajectories
        return trajectories


class StochKitODESolver(StochKitBaseSolver):
    name = "StochKitODESolver"
    """
    Abstract class for StochKit solver derived from the GillesPySolver class.
    This is generally used to set up the solver.

    Attributes
    ----------
    model : gillespy.Model
        The model on which the solver will operate.
    t : float
        The end time of the solver.
    number_of_trajectories : int
        The number of times to sample the chemical master equation. Each
        trajectory will be returned at the end of the simulation.
    increment : float
        The time step of the solution.
    seed : int
        The random seed for the simulation. Defaults to None.
    stochkit_home : str
        Path to stochkit. This is set automatically upon installation, but
        may be overwritten if desired.
    algorithm : str
        Already set to 'stochkit_ode.py'
    job_id : str
        If given, this will be the name of the solver run. Usually not set.
    debug : bool (False)
        Set to True to provide additional debug information about the
        simulation.
    """

    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, stochkit_home=None,
            algorithm='stochkit_ode.py',
            job_id=None, debug=False, profile=False, show_labels=False, **kwargs):
        return super().run(model, t, number_of_trajectories, increment, seed, stochkit_home,
                           algorithm, job_id, debug=debug, show_labels=show_labels)

    @classmethod
    def get_trajectories(cls, out_dir, debug=False, show_labels=False):
        if debug:
            print("StochKitODESolver.get_trajectories(out_dir={0}".format(out_dir))
        # Collect all the output data
        trajectories = []
        with open(os.path.join(out_dir, 'output.txt'), 'r') as fd:
            fd.readline()
            headers = fd.readline()
            fd.readline()
            data = [[float(x) for x in fd.readline().split()]]
            fd.readline()
            for line in fd:
                data.append([float(x) for x in line.split()])
        trajectories.append(np.array(data))
        if show_labels:
            return headers.split(), trajectories
        return trajectories
