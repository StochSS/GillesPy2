"""
An optional solver for simulating models which runs until desired convergence is met.

This solver serves primarily as a C wrapper to the ssa_c_solver.
It utilizes posix-threading to run multiple C_Solver processes, and stores
them in two different structures (evens and odds). The user designates
a desired "alpha" which correlates to a maximum distance determined by a
is calculated by a Kolmogorov-Smirnov distance comparison between
the evens and odds.

"""

import os #for getting directories for C++ files
import shutil #for deleting/copying files
import subprocess #For calling make and executing c solver
import inspect #for finding the Gillespy2 module path
import tempfile #for temporary directories
import numpy as np
import gillespy2
from .gillespySolver import GillesPySolver

GILLESPY_PATH = os.path.dirname(inspect.getfile(gillespy2))
GILLESPY_C_DIRECTORY = os.path.join(GILLESPY_PATH, 'c_base/')
GILLESPY_C_MULTI_DIRECTORY = os.path.join(GILLESPY_C_DIRECTORY, 'c_multi_solver')

def copy_files(destination):
    """
    Copy compile files to target destination
    Attributes
    ----------
    destination : str
        Target location for executable files.
    """

    src_files = os.listdir(GILLESPY_C_DIRECTORY)
    for src_file in src_files:
        src_file = os.path.join(GILLESPY_C_DIRECTORY, src_file)
        if os.path.isfile(src_file):
            shutil.copy(src_file, destination)
        else:
            if os.path.isdir(src_file):
                shutil.copytree(src_file, os.path.join(destination, 'c_multi_solver'))

def write_constants(outfile, model, reactions, species):
    """
    Write constants to simulation template.

    Attributes
    ----------
    outfile : str
        File created from template for compile.
    model : Model
        model to be simulated.
    reactions : dict
        model reactions for simulation
    species : dict
        model species for simulation
    """
    # Write mandatory constants
    outfile.write("const double vol = {};\n".format(model.volume))
    outfile.write("std :: string s_names[] = {")
    if species:
        # Write model species names.
        for i in range(len(species) - 1):
            outfile.write('"{}", '.format(species[i]))
        outfile.write('"{}"'.format(species[-1]))
        outfile.write("};\nuint populations[] = {")
        # Write initial populations.
        for i in range(len(species) - 1):
            outfile.write('{}, '.format(model.listOfSpecies[species[i]].initial_value))
        outfile.write('{}'.format(model.listOfSpecies[species[-1]].initial_value))
        outfile.write("};\n")
    if reactions:
        # Write reaction names
        outfile.write("std :: string r_names[] = {")
        for i in range(len(reactions) - 1):
            outfile.write('"{}", '.format(reactions[i]))
        outfile.write('"{}"'.format(reactions[-1]))
        outfile.write("};\n")
    for param in model.listOfParameters:
        outfile.write("const double {0} = {1};\n".format
                      (param, model.listOfParameters[param].value))


def write_propensity(outfile, model, reactions, species):
    """
    Write model propensities to simulation template
    Attributes
    ----------
    outfile : str
        File created from template for compile.
    model : Model
        model to be simulated.
    reactions : dict
        model reactions for simulation
    species : dict
        model species for simulation
    """
    for i in range(len(reactions)):
        propensity_function = model.listOfReactions[reactions[i]].propensity_function
        # Replace species references with array references
        for j in range(len(species)):
            propensity_function = propensity_function.replace(species[j], "state[{}]".format(j))
        # Write switch statement case for reaction
        outfile.write("""
        case {0}:
            return {1};
        """.format(i, propensity_function))


def write_reactions(outfile, model, reactions, species):
    """
    Write model reactions to simulation template.

    Attributes
    ----------
    outfile : str
       File created from template for compile.
    model : Model
        model to be simulated.
    reactions : dict
        model reactions for simulation
    species : dict
        model species for simulation
    """
    for i in range(len(reactions)):
        reaction = model.listOfReactions[reactions[i]]
        for j in range(len(species)):
            change = (reaction.products.get(model.listOfSpecies[species[j]], 0)) - (
                reaction.reactants.get(model.listOfSpecies[species[j]], 0))
            if change != 0:
                outfile.write("model.reactions[{0}].species_change[{1}] = {2};\n".format
                              (i, j, change))


def parse_output(results, number_timesteps, number_species):
    """
    Parses the results form the simulation executable into a
    3-dimensional numpy array with dimensions [trajectory][timestep][species]

    Attributes
    ----------
    results : str
        The simulation results returned by the executable as a string.
    number_timesteps : int
        Number of timesteps per trajectory
    number_species : int
        number of species present in the model
    """
    values = bytes(results).decode('utf-8')
    values = values.splitlines()
    number_of_trajectories = int(values[-1].rsplit()[-1])
    trajectory_base = np.empty((number_of_trajectories, number_timesteps, number_species+1))
    trajectory_n = 0
    timestep = 0
    lines_to_read = (number_of_trajectories * number_timesteps) + 1
    number_of_runs = values[-1]
    print(number_of_runs)
    for line in range(1, lines_to_read):
        if timestep == number_timesteps:
            trajectory_n += 1
            timestep = 0
        value = values[line].strip().split(" ")
        trajectory_base[trajectory_n, timestep, 0] = float(value[0])
        for species in range(number_species):
            trajectory_base[trajectory_n, timestep, species+1] = value[species+1]
        timestep += 1
    return trajectory_base

class SSACMultiSolver(GillesPySolver):
    """
        Optional Solver. Solves with basic SSA until desired convergence is met.

        Attributes
        ----------
        name : str
            The name by which this solver will be called.
        model : gillespy2.Model
            The model associated with this instance of the solver.
        output_directory : str
            desired output_directory for compiled C files
        delete_directory : bool
            if True, output_directory is removed after run
        number_of_processes : int
            number of processes to execute the simulations
        alpha : float
            desired maximum ks distance for convergence
        win_py_native : bool
            set to True if you are running this program on Windows with native Python
        """
    name = "SSACMultiSolver"
    """TODO"""

    def __init__(self, model=None, output_directory=None, delete_directory=True,
                 number_of_processes=4, alpha=0.01, win_py_native=False):
        super(SSACMultiSolver, self).__init__()
        self.compiled = False
        self.delete_directory = False
        self.model = model
        self.number_of_processes = number_of_processes
        self.alpha = alpha
        self.win_py_native = win_py_native
        self.simulation_data = None

        if alpha > 1 or self.alpha <= 0:
            raise gillespy2.InvalidAlphaError("Alpha must be between 0 and 1, 1 inclusive")

        if number_of_processes > 100:
            raise gillespy2.InvalidProcessesError("Number of Processes must be between 1 and 100")

        if self.model:
            # Create constant, ordered lists for reactions/species
            self.reactions = list(self.model.listOfReactions.keys())
            self.species = list(self.model.listOfSpecies.keys())

            if isinstance(output_directory, str):
                output_directory = os.path.abspath(output_directory)

                if isinstance(output_directory, str) and not os.path.isfile(output_directory):
                    self.output_directory = output_directory
                    self.delete_directory = delete_directory
                    if not os.path.isdir(output_directory):
                        # set up directory if needed
                        os.makedirs(self.output_directory)
            else:
                # Set up temporary directory
                self.temporary_directory = tempfile.TemporaryDirectory()
                self.output_directory = self.temporary_directory.name
            # copy files to directory
            copy_files(self.output_directory)
            # write template file
            self.write_template()
            # compile file
            self.compile()

    def __del__(self):
        if self.delete_directory and os.path.isdir(self.output_directory):
            shutil.rmtree(self.output_directory)

    def write_template(self, template_file='SimulationTemplate.cpp'):
        """
        Write simulation template file to be compiled.
        """
        # Open up template file for reading.
        with open(os.path.join(self.output_directory, template_file), 'r') as template:
            # Write simulation C++ file.
            template_keyword = "__DEFINE_"
            # Use same lists of model's species and reactions to maintain order
            with open(os.path.join(self.output_directory, 'UserSimulation.cpp'), 'w') as outfile:
                for line in template:
                    if line.startswith(template_keyword):
                        line = line[len(template_keyword):]
                        if line.startswith("CONSTANTS"):
                            write_constants(outfile, self.model, self.reactions, self.species)
                        if line.startswith("PROPENSITY"):
                            write_propensity(outfile, self.model, self.reactions, self.species)
                        if line.startswith("REACTIONS"):
                            write_reactions(outfile, self.model, self.reactions, self.species)
                    else:
                        outfile.write(line)

    def compile(self):
        """
        Compile the simulation files into Linux executables
        """
        multi_path = os.path.join(self.output_directory, 'c_multi_solver')
        bash_path = os.path.join(os.environ['SystemRoot'], 'System32', 'bash.exe')
        bash_path2 = os.path.join(os.environ['SystemRoot'], 'SysNative', 'bash.exe')
        if self.win_py_native:
            bash_path = bash_path2
        # Use makefile.
        comp_arg = '{0} -c make'.format(bash_path)
        subprocess.run(
            ["make", "-C", self.output_directory, 'cleanSimulation'], stdout=subprocess.PIPE)
        built = subprocess.run(
            ["make", "-C", self.output_directory, 'UserSimulation'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if os.name == 'nt':
            multi_path = multi_path.replace('\\', '/')
            built_multi = subprocess.run(
                comp_arg, cwd=multi_path, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, shell=False)
        else:
            built_multi = subprocess.run(
                ["make", "-C", multi_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #print(built_multi.stdout)
        # Use makefile.
        if built.returncode == 0 and built_multi.returncode == 0:
            self.compiled = True
        else:
            print("Error encountered while compiling file:\n"
                  "built Return code: {0}.\nError:\n{1}\n"
                  .format(built.returncode, built.stderr))
            print("Error encountered while compiling file:\n"
                  "built_multi Return code: {0}.\nError:\n{1}\n"
                  .format(built_multi.returncode, built_multi.stderr))

    def run(self=None, model=None, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False, stochkit_home=None):
        """
            Runs the solver.

            Attributes
            ----------
            model : Model
                model to perform the solver on
            t : int
                total run time
            increment : float
                time step size for along time span t
            seed : int
                The random seed for the simulation. Optional, defaults to None.
            debug : bool (False)
                Set to True to provide additional debug information about the
                simulation.
            show_labels : bool (True)
                Use names of species as index of result object rather than position numbers.
            """

        if self is None:
            self = SSACMultiSolver(model)

        if self.compiled:
            number_timesteps = int(t // increment)
            if t % increment > 0:
                number_timesteps += 1
            # Execute simulation.
            multi_path = os.path.join(self.output_directory, 'c_multi_solver/')
            multi_path = multi_path.replace('\\', '/')
            args = self.getargs(seed, number_timesteps, t)
            simulation = subprocess.run(
                args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=multi_path, shell=False)
            # Parse/return results.

            if simulation.returncode == 0:
                trajectory_base = parse_output(
                    simulation.stdout, number_timesteps, len(self.species))
                # Format results
                if show_labels:
                    self.simulation_data = []
                    for trajectory in range(len(trajectory_base)):
                        data = {}
                        data['time'] = trajectory_base[trajectory, :, 0]
                        for i in range(len(self.species)):
                            data[self.species[i]] = trajectory_base[trajectory, :, i+1]
                        self.simulation_data.append(data)
                else:
                    self.simulation_data = trajectory_base
            else:
                print("Error encountered while running simulation C++ file:\n"
                      "Return code: {0}.\nError:\n{1}\n"
                      .format(simulation.returncode, simulation.stderr))

        return self.simulation_data

    def getargs(self, seed, number_timesteps, t):
        """
        Assign correct subprocess calls for current system and environment
        """
        bash_path = os.path.join(os.environ['SystemRoot'], 'System32', 'bash.exe')
        bash_path2 = os.path.join(os.environ['SystemRoot'], 'SysNative', 'bash.exe')

        if self.win_py_native:
            bash_path = bash_path2


        multi_path = os.path.join(self.output_directory, 'c_multi_solver/')
        multi_path = multi_path.replace('\\', '/')
        if os.name == 'nt':
            args = '{0} -c "./c_multi_solver ../UserSimulation.exe {1} {2} {3} {4} {5}"'.format(
                bash_path, str(self.number_of_processes), str(len(self.species)),
                str(number_timesteps), str(t), str(self.alpha))
            if isinstance(seed, int):
                args.append(str(seed))

        else:
            args = "./c_multi_solver ../UserSimulation {1} {2} {3} {4} {5}".format(
                str(self.number_of_processes), str(len(self.species)),
                str(number_timesteps), str(t)), str(self.alpha)
            if isinstance(seed, int):
                args.append('-seed')
                args.append(str(seed))
        return args
