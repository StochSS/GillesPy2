import gillespy2
from .gillespySolver import GillesPySolver

def write_constants(outfile, model, t, number_of_trajectories, increment, seed, reactions, species):
    number_timesteps = int(t/increment)
    #Write mandatory constants
    outfile.write("""
const uint number_trajectories = {0};
const uint number_timesteps = {1};
const double end_time = {2};
const double vol = {3};
int random_seed;
""".format(number_of_trajectories, number_timesteps, end_time, model.volume))
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
        for j in range(len(species)):
            propensity_function = propensity_function.replace(species[j], "state[{}]".format(j))
        outfile.write("""

        case {0}:
            return {1};

        """.format(i, propensity_function))
        
class SSACSolver(GillesPySolver):
    """TODO"""
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        #Open up template file for reading.
        with open('SimulationTemplate.cpp', 'r') as template:
            #Write simulation C++ file.
            template_keword = "__DEFINE_"
            #Use same lists of model's species and reactions to maintain order
            reactions = model.listOfReactions.keys()
            species = model.listOfSpecies.keys()
            with open('UserSimulation.cpp', 'w') as outfile:
                for line in template:
                    if line.startswith(template_keyword):
                        line = line[len(template_keyword):]
                        if line.startswith("CONSTANTS"):
                            write_constants(outfile, model, t, number_of_trajectories, increment, seed)
                        if line.startswith("PROPENSITY"):
                            write_propensity(outfile, model)
                    else:
                        outfile.write(line)
                #Write propensity function.
                pass
        #Use makefile.
        #Execute simulation.
        #Parse/return results.
        pass
