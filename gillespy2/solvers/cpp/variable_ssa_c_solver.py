import gillespy2
from gillespy2.core import Model, Reaction, gillespyError, GillesPySolver, log
from gillespy2.solvers.utilities.utilities import species_parse
import signal, time #for solver timeout implementation
import os #for getting directories for C++ files
import shutil #for deleting/copying files
import subprocess #For calling make and executing c solver
import inspect #for finding the Gillespy2 module path
import tempfile #for temporary directories
import numpy as np

GILLESPY_PATH = os.path.dirname(inspect.getfile(gillespy2))
GILLESPY_C_DIRECTORY = os.path.join(GILLESPY_PATH, 'solvers/cpp/c_base')


def find_time(array,value):
    """
    Finds the index of the closest value in the array parameter, to the value parameter
    :param array: results['time'] array, input to find index of closest 'value' parameter
    :type array: numpy.ndarray
    :param value: Value in which to find the closest values index in the array parameter.
    :type value: float
    :return: Integer index, the index of the closest value to 'value' parameter.
    """
    index = np.searchsorted(array, value, side="left")
    return index

def _copy_files(destination):
    src_files = os.listdir(GILLESPY_C_DIRECTORY)
    for src_file in src_files:
        src_file = os.path.join(GILLESPY_C_DIRECTORY, src_file)
        if os.path.isfile(src_file):
            shutil.copy(src_file, destination)


def _write_variables(outfile, model, reactions, species, parameters, parameter_mappings, resume=None):
    outfile.write("double V = {};\n".format(model.volume))
    outfile.write("std :: string s_names[] = {")
    if len(species) > 0:
        #Write model species names.
        for i in range(len(species)-1):
            outfile.write('"{}", '.format(species[i]))
        outfile.write('"{}"'.format(species[-1]))
        outfile.write("};\nunsigned int populations[] = {")
        #Write initial populations.
        for i in range(len(species) - 1):
            # If resuming
            if not (resume is None):
                outfile.write('{}, '.format(int(resume[species[i]][-1])))
            else:
                outfile.write('{}, '.format(int(model.listOfSpecies[species[i]].initial_value)))
        if not (resume is None):
            if isinstance(resume, np.ndarray):
                outfile.write('{}'.format(int(resume[0][-1][-1])))
            else:
                outfile.write('{}'.format(int(resume[species[-1]][-1])))
        else:
            outfile.write('{}'.format(int(model.listOfSpecies[species[-1]].initial_value)))
        outfile.write("};\n")
    if len(reactions) > 0:
        #Write reaction names
        outfile.write("std :: string r_names[] = {")
        for i in range(len(reactions)-1):
            outfile.write('"{}", '.format(reactions[i]))
        outfile.write('"{}"'.format(reactions[-1]))
        outfile.write("};\n")
    for param in parameters:
        if param != 'vol':
            outfile.write("double {0} = {1};\n".format(parameter_mappings[param], model.listOfParameters[param].value))

def _update_parameters(outfile, model, parameters, parameter_mappings):
    for param in parameters:
        if param != 'vol':
            outfile.write('       arg_stream >> {};\n'.format(parameter_mappings[param]))
        else:
            outfile.write('       arg_stream >> V;\n')

def _write_propensity(outfile, model, species_mappings, parameter_mappings, reactions):
    for i in range(len(reactions)):
        # Write switch statement case for reaction
        outfile.write("""
        case {0}:
            return {1};
        """.format(i, model.listOfReactions[reactions[i]].sanitized_propensity_function(species_mappings, parameter_mappings)))


def _write_reactions(outfile, model, reactions, species):
    customrxns = {}
    for i in range(len(reactions)):
        reaction = model.listOfReactions[reactions[i]]
        if reaction.type == 'customized':
            customrxns[i] = species_parse(model, reaction.propensity_function)
        for j in range(len(species)):
            change = (reaction.products.get(model.listOfSpecies[species[j]], 0)) - (
                reaction.reactants.get(model.listOfSpecies[species[j]], 0))
            if change != 0:
                outfile.write("model.reactions[{0}].species_change[{1}] = {2};\n".format(i, j, change))

    for i in customrxns.keys():
        for j in range(len(reactions)):
            if i == j:
                continue
            if any(elem in customrxns[i] for elem in list(model.listOfReactions[reactions[j]].reactants)) or \
                    any(elem in customrxns[i] for elem in list(model.listOfReactions[reactions[j]].products)):
                outfile.write("model.reactions[{0}].affected_reactions.push_back({1});\n".format(i, j))

def _parse_binary_output(results_buffer, number_of_trajectories, number_timesteps, number_species,pause=False):
    trajectory_base = np.empty((number_of_trajectories, number_timesteps, number_species+1))
    step_size = number_species * number_of_trajectories + 1 #1 for timestep
    data = np.frombuffer(results_buffer, dtype=np.float64)
    #Timestopped is added to the end of the data, when a simulation completes or is paused
    if pause:
        timeStopped = data[-1]
    else:
        timeStopped = 0
    assert(len(data) == (number_of_trajectories*number_timesteps*number_species + number_timesteps)+1)
    for timestep in range(number_timesteps):
        index = step_size * timestep
        trajectory_base[:, timestep, 0] = data[index]
        index += 1
        for trajectory in range(number_of_trajectories):
            for species in range(number_species):
                trajectory_base[trajectory, timestep, 1 + species] = data[index + species]
            index += number_species
    return trajectory_base, timeStopped


class VariableSSACSolver(GillesPySolver):
    name = "VariableSSACSolver"
    def __init__(self, model=None, output_directory=None, delete_directory=True, resume=None):
        super(VariableSSACSolver, self).__init__()
        self.__compiled = False
        self.delete_directory = False
        self.model = model
        self.resume = resume
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
                raise gillespyError.DirectoryError("Errors encountered while setting up directory for Solver C++ files.")
            _copy_files(self.output_directory)
            self.__write_template()
            self.__compile()
        
    def __del__(self):
        if self.delete_directory and os.path.isdir(self.output_directory):
            shutil.rmtree(self.output_directory)
        
    def __write_template(self, template_file='VariableSimulationTemplate.cpp'):
        # Open up template file for reading.
        with open(os.path.join(self.output_directory, template_file), 'r') as template:
            # Write simulation C++ file.
            template_keyword = "__DEFINE_"
            # Use same lists of model's species and reactions to maintain order
            with open(os.path.join(self.output_directory, 'UserSimulation.cpp'), 'w') as outfile:
                for line in template:
                    if line.startswith(template_keyword):
                        line = line[len(template_keyword):]
                        if line.startswith("VARIABLES"):
                            _write_variables(outfile, self.model, self.reactions, self.species, self.parameters,
                                             self.parameter_mappings,self.resume)
                        if line.startswith("PROPENSITY"):
                            _write_propensity(outfile, self.model, self.species_mappings, self.parameter_mappings,
                                              self.reactions)
                        if line.startswith("REACTIONS"):
                            _write_reactions(outfile, self.model, self.reactions, self.species)
                        if line.startswith("PARAMETER_UPDATES"):
                            _update_parameters(outfile, self.model, self.parameters, self.parameter_mappings)
                    else:
                        outfile.write(line)

    def __compile(self):
        # Use makefile.
        if self.resume:
            if self.resume[0].model != self.model:
                raise gillespyError.ModelError('When resuming, one must not alter the model being resumed.')
            else:
                built = subprocess.run(["make", "-C", self.output_directory, 'UserSimulation'], stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        else:
            try:
                cleaned = subprocess.run(["make", "-C", self.output_directory, 'cleanSimulation'],
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                built = subprocess.run(["make", "-C", self.output_directory, 'UserSimulation'],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            except KeyboardInterrupt:
                log.warning("Solver has been interrupted during compile time, unexpected behavior may occur.")

        if built.returncode == 0:
            self.__compiled = True
        else:
            raise gillespyError.BuildError("Error encountered while compiling file:\nReturn code: {0}."
                                           "\nError:\n{1}\n{2}\n".format(built.returncode, built.stdout.decode
            ('utf-8'),built.stderr.decode('utf-8')))

    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'timeout', 'increment', 'seed', 'debug', 'profile', 'variables')

    def run(self=None, model=None, t=20, number_of_trajectories=1, timeout=0,
            increment=0.05, seed=None, debug=False, profile=False, variables={}, resume=None, **kwargs):
        pause = False
        if resume is not None:
            if t < resume['time'][-1]:
                raise gillespyError.ExecutionError(
                    "'t' must be greater than previous simulations end time, or set in the run() method as the "
                    "simulations next end time")

        if resume is not None:
            self = VariableSSACSolver(model, resume=resume)
        else:
            if self is None or self.model is None:
                self = VariableSSACSolver(model)

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
                'Could not run Model.  SBML Feature: {} not supported by SSACSolver.'.format(detected_features))

        if not isinstance(variables, dict):
            raise gillespyError.SimulationError(
                'argument to variables must be a dictionary.')
        for v in variables.keys():
            if v not in self.species+self.parameters:
                raise gillespyError.SimulationError('Argument to variable "{}" \
                is not a valid variable.  Variables must be model species or parameters.'.format(v))
                
        if self.__compiled:
            populations = ''
            parameter_values = ''
            # Update Species Initial Values
            for i in range(len(self.species)-1):
                if self.species[i] in variables:
                    populations += '{} '.format(int(variables[self.species[i]]))
                else:
                    if resume is not None:
                        populations += '{} '.format(int(resume[self.species[i]][-1]))
                    else:
                        populations += '{} '.format(int(model.listOfSpecies[self.species[i]].initial_value))
            if self.species[-1] in variables:
                populations += '{}'.format(int(variables[self.species[-1]]))
            else:
                if resume is not None:
                    populations += '{} '.format(int(resume[self.species[-1]][-1]))
                else:
                    populations += '{}'.format(int(model.listOfSpecies[self.species[-1]].initial_value))
            # Update Parameter Values
            for i in range(len(self.parameters)-1):
                if self.parameters[i] in variables:
                    parameter_values += '{} '.format(variables[self.parameters[i]])
                else:
                    if self.parameters[i] == 'vol':
                        parameter_values +='{} '.format(model.volume)
                    else:
                        parameter_values +='{} '.format(model.listOfParameters[self.parameters[i]].expression)
            if self.parameters[-1] in variables:
                parameter_values += '{}'.format(variables[self.parameters[-1]])
            else:
                if self.parameters[i] == 'vol':
                    parameter_values += '{}'.format(model.volume)
                else:
                    parameter_values += '{}'.format(model.listOfParameters[self.parameters[-1]].expression)
            self.simulation_data = None

            if resume is not None:
                t = abs(t - int(resume['time'][-1]))

            number_timesteps = int(round(t/increment + 1))
            # Execute simulation.
            args = [os.path.join(self.output_directory, 'UserSimulation'), 
                    '-trajectories', str(number_of_trajectories), 
                    '-timesteps', str(number_timesteps), 
                    '-end', str(t),
                    '-initial_values', populations,
                    '-parameters', parameter_values]
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
                        raise ModelError("seed must be a positive integer")

            #begin subprocess c simulation with timeout (default timeout=0 will not timeout)
            with subprocess.Popen(args, stdout=subprocess.PIPE, start_new_session=True) as simulation:
                try:
                    if timeout > 0:
                        stdout, stderr = simulation.communicate(timeout=timeout)
                    else:
                        stdout, stderr = simulation.communicate()
                    return_code = simulation.wait()
                except KeyboardInterrupt:
                    os.killpg(simulation.pid, signal.SIGINT)  # send signal to the process group
                    stdout, stderr = simulation.communicate()
                    pause = True
                    return_code = 33
                except subprocess.TimeoutExpired:
                        os.killpg(simulation.pid, signal.SIGINT) #send signal to the process group
                        stdout, stderr = simulation.communicate()
                        pause = True
                        return_code = 33
 
            # Parse/return results.
            if return_code in [0, 33]:
                trajectory_base, timeStopped = _parse_binary_output(stdout, number_of_trajectories, number_timesteps,
                                                                    len(self.species),pause=pause)
                if self.model.tspan[2]-self.model.tspan[1] == 1:
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
            # If simulation was paused/KeyboardInterrupt
            if timeStopped != 0:
                cutoff = find_time(self.simulation_data[0]['time'], timeStopped)
                if cutoff == 0 or cutoff == 1:
                    log.warning('You have paused the simulation too early, and no points have been calculated past'
                                ' initial values. A graphic display will not produce expected results.')
                else:
                    cutoff -= 1
                for i in self.simulation_data[0]:
                    self.simulation_data[0][i] = self.simulation_data[0][i][:cutoff]

            if resume is not None:
                resumeTime = float(resume['time'][-1])
                step = resumeTime - resume['time'][-2]
                if timeStopped == 0:
                    timeSpan = np.arange(resumeTime, t + resumeTime + step, step)
                else:
                    timeSpan = np.arange(resumeTime + step, timeStopped + resumeTime + step, step)
                self.simulation_data[0]['time'] = timeSpan

            if resume is not None:
                # If resuming, combine old pause with new data, and delete any excess null data
                for i in self.simulation_data[0]:
                    oldData = resume[i]
                    newData = self.simulation_data[0][i]
                    self.simulation_data[0][i] = np.concatenate((oldData, newData), axis=None)
                if len(self.simulation_data[0]['time']) != len(self.simulation_data[0][i]):
                    self.simulation_data[0]['time'] = self.simulation_data[0]['time'][:-1]

        return self.simulation_data, return_code

