import gillespy2
from .gillespySolver import GillesPySolver

def write_constants(outfile, model, t, number_of_trajectories, increment, seed):
    number_timesteps = int(t/increment)
    outfile.write("""
const uint number_trajectories = {0};
const uint number_timesteps = {1};
const double end_time = {2};
const double volume = {3};
int random_seed;
""".format(number_of_trajectories, number_timesteps, end_time, model.volume))
    #Write seed
    if isinstance(seed, int):
        outfile.write("random_seed = {};\nseed_time = false;\n".format(seed))        
    #    Write model species names.
    outfile.write("std :: string s_names[] = {");
    for species in model.listOfSpecies:
        outfile.write('"{}"'.format(species))
    #    Write initial populations.
         


class SSACSolver(GillesPySolver):
    """TODO"""
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        #Open up template file for reading.
        with open('SimulationTemplate.cpp', 'r') as template:
            #Write simulation C++ file.
            template_keword = "__DEFINE_"
            with open('UserSimulation.cpp', 'w') as outfile:
                for line in template:
                    if line.startswith(template_keyword):
                        line = line[len(template_keyword):]
                        if line.startswith("CONSTANTS"):
                            write_constants(outfile, model, t, number_of_trajectories, increment, seed)
                    else:
                        outfile.write(line)
                #Write propensity function.
                pass
        #Use makefile.
        #Execute simulation.
        #Parse/return results.
        pass
