import gillespy2
from .gillespySolver import GillesPySolver
import os #for getting directories for C++ filesy
import subprocess #For calling make and executing c solver
import inspect #for finding the Gillespy2 module path
import numpy as np
import math

def write_constants(outfile, model, t, number_of_trajectories, number_timesteps, seed, reactions, species):
    #Write mandatory constants
    outfile.write("""
const uint number_trajectories = {0};
const uint number_timesteps = {1};
const double end_time = {2};
const double vol = {3};
int random_seed;
""".format(number_of_trajectories, number_timesteps, t, model.volume))
    #Write seed
    if isinstance(seed, int):
        outfile.write("random_seed = {};\nseed_time = false;\n".format(seed))        
    outfile.write("std :: string s_names[] = {");
    if len(species) > 0:
        #Write model species names.
        for i in range(len(species)-1):
            outfile.write('"{}", '.format(species[i]))
        outfile.write('"{}"'.format(species[-1]))
        outfile.write("};\nuint populations[] = {")
        #Write initial populations.
        for i in range(len(species)-1):
            outfile.write('{}, '.format(model.listOfSpecies[species[i]].initial_value))
        outfile.write('{}'.format(model.listOfSpecies[species[-1]].initial_value))
        outfile.write("};\n")
    if len(reactions) > 0:
        #Write reaction names
        outfile.write("std :: string r_names[] = {")
        for i in range(len(reactions)-1):
            outfile.write('"{}", '.format(reactions[i]))
        outfile.write('"{}"'.format(reactions[-1]))
        outfile.write("};\n")
    for param in model.listOfParameters:
        outfile.write("const double {0} = {1};\n".format(param, model.listOfParameters[param].value))


def write_propensity(outfile, model, reactions, species):
    for i in range(len(reactions)):
        propensity_function = model.listOfReactions[reactions[i]].propensity_function
        #Replace species references with array references
        for j in range(len(species)):
            propensity_function = propensity_function.replace(species[j], "state[{}]".format(j))
        #Write switch statement case for reaction
        outfile.write("""

        case {0}:
            return {1};

        """.format(i, propensity_function))


def write_reactions(outfile, model, reactions, species):
    for i in range(len(reactions)):
        reaction = model.listOfReactions[reactions[i]]
        for j in range(len(species)):
            change = (reaction.products.get(model.listOfSpecies[species[j]], 0)) - (reaction.reactants.get(model.listOfSpecies[species[j]], 0))
            if change != 0:
                outfile.write("model.reactions[{0}].species_change[{1}] = {2};\n".format(i, j, change))


def parse_output(results, number_of_trajectories, number_timesteps, number_species):
    trajectory_base = np.empty((number_of_trajectories, number_timesteps, number_species+1))
    for timestep in range(number_timesteps):
        values = results[timestep].split(" ")
        time = float(values[0])
        index = 1
        for trajectory in range(number_of_trajectories):
            trajectory_base[trajectory, timestep, 0] = time
            for species in range(number_species):
                trajectory_base[trajectory, timestep, 1 + species] = float(values[index+species])
            index += number_species
    return trajectory_base

class SSACSolver(GillesPySolver):
    """TODO"""
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        GILLESPY_PATH = os.path.dirname(inspect.getfile(gillespy2))
        GILLESPY_C_DIRECTORY = os.path.join(GILLESPY_PATH, 'c_base/')
        self.simulation_data = None
        number_timesteps = int(math.ceil(t/increment))
        #Open up template file for reading.
        with open(os.path.join(GILLESPY_C_DIRECTORY,'SimulationTemplate.cpp'), 'r') as template:
            #Write simulation C++ file.
            template_keyword = "__DEFINE_"
            #Use same lists of model's species and reactions to maintain order
            reactions = list(model.listOfReactions.keys())
            species = list(model.listOfSpecies.keys())
            with open(os.path.join(GILLESPY_C_DIRECTORY, 'UserSimulation.cpp'), 'w') as outfile:
                for line in template:
                    if line.startswith(template_keyword):
                        line = line[len(template_keyword):]
                        if line.startswith("CONSTANTS"):
                            write_constants(outfile, model, t, number_of_trajectories, number_timesteps, seed, reactions, species)
                        if line.startswith("PROPENSITY"):
                            write_propensity(outfile, model, reactions, species)
                        if line.startswith("REACTIONS"):
                            write_reactions(outfile, model, reactions, species)
                    else:
                        outfile.write(line)
        #Use makefile.
        cleaned = subprocess.run(["make", "-C", GILLESPY_C_DIRECTORY, 'cleanSimulation'], stdout=subprocess.PIPE)
        built = subprocess.run(["make", "-C", GILLESPY_C_DIRECTORY, 'UserSimulation'], stdout=subprocess.PIPE)
        if built.returncode == 0:
            #Execute simulation.
            simulation = subprocess.run([os.path.join(GILLESPY_C_DIRECTORY, 'UserSimulation')], stdout=subprocess.PIPE)
            #Parse/return results.
            if simulation.returncode == 0:
                results = simulation.stdout.decode('utf-8').split('\n')
                trajectory_base = parse_output(results, number_of_trajectories, number_timesteps, len(species))
                #Format results
                if show_labels:
                    self.simulation_data = []
                    for trajectory in range(number_of_trajectories):
                        data = {}
                        data['time'] = trajectory_base[trajectory,:,0]
                        for i in range(len(species)):
                            data[species[i]] = trajectory_base[trajectory, :, i]
                        self.simulation_data.append(data)
                else:
                    self.simulation_data = trajectory_base
        return self.simulation_data

