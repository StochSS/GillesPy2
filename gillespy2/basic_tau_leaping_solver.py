import gillespy2
from gillespySolver import GillesPySolver
import random
import math
import numpy

class BasicTauLeapingSolver(GillesPySolver):
    """ TODO
    """
    
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        self.simulation_data = []
 
        # initialize everything
        curr_state = {}
        propensity = {}
        results = {}
        poisson_values = {}
        expected_changes = {}

        for traj_num in range(number_of_trajectories):
            for s in model.listOfSpecies:   #Initialize Species population
                curr_state[s] = model.listOfSpecies[s].initial_value    
                results[s]=[]

            curr_state['vol'] = model.volume
            results['time'] = []
            curr_time = 0

            for p in model.listOfParameters:
                curr_state[p] = model.listOfParameters[p].value        
        
            # run the loop 
            while(curr_time < t):
                # calculate propensities
                for r in model.listOfReactions: 
                    propensity[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
                    prop_sum += propensity[r]    

                # pick tau.
                # fixed tau for nau
                tau = 0.2
                curr_time += tau

                # generate Poisson distributions
                for r in model.listOfReactions:
                    poisson_values[r] = np.random.poisson(propensity[r],tau)

                # update curr_state based on the expected changes of each reaction
                for r in model.listOfReactions:
                    for react in model.listOfReactions[r].reactants:
                        curr_state[react] -= poisson_values[r]*model.listOfReactions[r].reactants[react]
                    for prod in model.listOfReactions[r].products:
                        curr_state[prod] += poisson_values[r].model.listOfReactions[r].products[prod]

                # append curr_state to results
                results['time'].append(curr_time)
                for s in model.listOfSpecies:
                    results[s].append(curr_state[s])

            return results

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if show_labels:
            return self.simulation_data
       # else:
            #TODO: need to account for 'show_labels'




