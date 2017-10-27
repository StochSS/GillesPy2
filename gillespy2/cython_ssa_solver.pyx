import gillespy2
from gillespySolver import GillesPySolver
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

cdef struct CythonSpecies:
    int index
    char* name
    double initial_population

cdef struct CythonReaction:
    int index
    char* name
    double propensity
    int *affected_reactions


def CythonSSASolver(GillesPySolver):
    """ TODO
    """

    
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
    #convert dictionary of species to species array
        cdef int number_species = len(list(model.listOfSpecies.keys()))
        cdef CythonSpecies *species = <CythonSpecies*> malloc(number_species * sizeof(CythonSpecies)) 
        cdef int i = 0
        for name, spec in model.listOfSpecies.items():
            species[i].index = i
            species[i].name = name
            species[i].initial_population = spec.initial_value
            i += 1
    #convert dictionary of reactions to reactions array
        cdef int number_reactions = len(list(model.listOfReactions.keys()))
        cdef CythonReaction *reactions = <CythonReaction*> malloc(number_reactions * sizeof(CythonReaction))
        i = 0
        for name, reaction in model.listOfReactions.items():
            reactions[i].index = i
            reactions[i].name = name
            reactions[i].affected_reactions = <int*> malloc(number_reactions*sizeof(int))
            i += 1
    #convert propensity functions now
        #create dictionary of all constant parameters for propensity evaluation
        parameters = {'vol' : model.volume}
        for paramName, param in model.listOfParameters.items():
            parameters[paramName] = param.value
        propensity_functions = [r.propensity_function for r in model.listOfReactions.values()]
        #create an array mapping reactions to species modified
        propensity_functions = []
        cdef np.ndarray[np.float64_t, ndim=2] species_changes = np.zeros([number_reactions, number_species])
        #pre-evaluate propensity equations from strings:
        for i in range(number_reactions):
            reaction = model.listOfReactions[reactions[i].name]
            #replace all references to species with array indices
            for j in range(number_species):
                spec = model.listOfSpecies[species[j].name]
                species_changes[i][j] = reaction.products.get(spec,0) - reaction.reactants.get(spec, 0)
                propensity_functions[i] = propensity_functions[i].replace(species[j].name, 'x[{0}]'.format(j))
            propensity_functions[i] = eval('lambda x:'+propensity_functions[i], parameters)
        #clean up
        for i in range(number_reactions):
            free(reactions[i].affected_reactions)
            
        free(reactions)
        free(species)
        
