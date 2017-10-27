import gillespy2
from gillespySolver import GillesPySolver
cimport numpy as np
from libc.stdlib cimport malloc, free

cdef struct CythonSpecies:
    int index
    char* name
    double initial_population

cdef struct CythonReaction:
    int index
    char* name
    double propensity
    double *propensity_function(np.ndarray[np.float64_t, ndim=1] state)
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
        #clean up
        for i in range(number_reactions):
            free(reactions[i].affected_reactions)
            
        free(reactions)
        free(species)
        
