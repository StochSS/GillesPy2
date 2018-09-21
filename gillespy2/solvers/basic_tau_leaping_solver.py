import gillespy2
from .gillespySolver import GillesPySolver
import random
import math
import numpy as np

class BasicTauSolver(GillesPySolver):
    name = "BasicTauSolver"
    
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        self.simulation_data = []
 
        curr_state = {}
        propensities = {}
        results = {}
        poissonValues = {}

        for traj_num in range(number_of_trajectories):
            for species in model.listOfSpecies:   #Initialize Species population
                curr_state[species] = model.listOfSpecies[species].initial_value    
                results[species]=[]

            curr_state['vol'] = model.volume
            results['time'] = []
            currentTime = 0
            outputTime = 0      
            nextTime = 0

            for parameter in model.listOfParameters:
                curr_state[parameter] = model.listOfParameters[parameter].value        
       
            # run the algorithm
            while(currentTime < t):
                while (currentTime >= outputTime) and (outputTime < t):
                    results['time'].append(outputTime)
                    for species in model.listOfSpecies:
                        results[species].append(curr_state[species])
                    outputTime += increment

                # evaluate propensities
                for reaction in model.listOfReactions: 
                    propensities[reaction] = eval(model.listOfReactions[reaction].propensity_function, curr_state)
   
                # select tau
                tau = 0.05
                # don't leap past next output time
                if currentTime + tau > outputTime:
                    tau = outputTime - currentTime
                    nextTime = outputTime
                else:
                    nextTime = currentTime + tau

                # calculate poisson distribution
                for reaction in model.listOfReactions:
                    poissonValues[reaction] = np.random.poisson(propensities[reaction]*tau) 

                # append changes to curr_state
                for reaction in model.listOfReactions:
                    for reactant in model.listOfReactions[reaction].reactants:
                        curr_state[str(reactant)] -= model.listOfReactions[reaction].reactants[reactant]*poissonValues[reaction]
                    for product in model.listOfReactions[reaction].products:
                        curr_state[str(product)] += model.listOfReactions[reaction].products[product]*poissonValues[reaction]

                # update the time
                currentTime = nextTime

            return results

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if show_labels:
            return self.simulation_data
       # else:
            #TODO: need to account for 'show_labels'




