import gillespy2
from gillespySolver import GillesPySolver
import numpy as np
import random
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.math cimport log

cdef struct CythonSpecies:
    int index
    char* name

cdef struct CythonReaction:
    int index
    double propensity
    int *affected_reactions


class CythonSSASolver(GillesPySolver):
    """ TODO
    """

    
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        self.simulation_data = []
        #convert dictionary of species to species array
        species = list(model.listOfSpecies.keys())
        cdef int number_species = len(species)
        #set timespan for simulation(s)
        timeline = np.linspace(0,t, (t//increment+1))
        #allocate memory for trajectories
        cdef np.ndarray[np.float64_t, ndim=3] trajectories = np.zeros((number_of_trajectories, timeline.size, number_species + 1))
        trajectories[:,:,0] = timeline
        cdef int i = 0, j
        for i in range(number_species):
            trajectories[:,:,i+1] = model.listOfSpecies[species[i]].initial_value
    #convert dictionary of reactions to reactions array
        cdef int number_reactions = len(list(model.listOfReactions.keys()))
        cdef CythonReaction *reactions = <CythonReaction*> malloc(number_reactions * sizeof(CythonReaction))
        i = 0
        reaction_names = list(model.listOfReactions.keys())
        for name in model.listOfReactions:
            reactions[i].index = i
            reactions[i].affected_reactions = <int*> malloc(number_reactions*sizeof(int))
            i += 1
    #convert propensity functions now
        #create dictionary of all constant parameters for propensity evaluation
        parameters = {'vol' : model.volume}
        for paramName, param in model.listOfParameters.items():
            parameters[paramName] = param.value
        propensity_functions = [r.propensity_function for r in model.listOfReactions.values()]
        cdef np.ndarray[np.float64_t, ndim=2] species_changes = np.zeros((number_reactions, number_species))
        #pre-evaluate propensity equations from strings:
        for i in range(number_reactions):
            reaction = model.listOfReactions[reaction_names[i]]
            #replace all references to species with array indices
            for j in range(number_species):
                spec = model.listOfSpecies[species[j]]
                species_changes[i][j] = reaction.products.get(spec,0) - reaction.reactants.get(spec, 0)
                propensity_functions[i] = propensity_functions[i].replace(species[j], 'x[{0}]'.format(j))
            propensity_functions[i] = eval('lambda x:'+propensity_functions[i], parameters)
        #begin simulating each trajectory
        cdef double current_time, propensity_sum, cumulative_sum
        cdef np.ndarray[np.float64_t, ndim=1] current_state = np.zeros((number_species))
        for i in range(number_of_trajectories):
            current_time = 0
            number_entries = 0
            np.copyto(current_state, trajectories[i,0,1:])
            for j in range(number_reactions):
                reactions[j].propensity = propensity_functions[j](current_state)

            while number_entries < timeline.size:
                propensity_sum = 0
                for j in range(number_reactions):
                    propensity_sum += reactions[j].propensity
                if propensity_sum <= 0:
                    trajectories[i,number_entries:,1:] = current_state
                    break
                cumulative_sum = random.random() * propensity_sum
                current_time -= log(random.random()) / propensity_sum
                while number_entries < timeline.size and timeline[number_entries] <= current_time:
                    trajectories[i,number_entries, 1:] = current_state
                    number_entries += 1
                for j in range(number_reactions):
                    cumulative_sum -= reactions[j].propensity
                    if cumulative_sum <= 0:
                        current_state += species_changes[j]
                        for j in range(number_reactions):
                            reactions[j].propensity = propensity_functions[j](current_state)
                        break
            #assemble complete simulation data in format specified
            if show_labels:
                data = {'time' : timeline}
                for j in range(number_species):
                    data[species[j]] = trajectories[i,:,j+1]
                self.simulation_data.append(data)
            else:
                self.simulation_data.append(trajectories[i])
        #clean up
        for i in range(number_reactions):
            free(reactions[i].affected_reactions)
        free(reactions)
        return self.simulation_data
        
