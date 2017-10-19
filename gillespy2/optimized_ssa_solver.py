import gillespy2
from gillespySolver import GillesPySolver
import random
import math
import numpy as np

import numba

class SSASolver(GillesPySolver):
    """ TODO
    """

    def format_trajectories(simulation_data):
        out_data = []
        sorted_keys = sorted(simulation_data[0])
        sorted_keys.remove('time')
        for trajectory in simulation_data:
            columns = [np.vstack((trajectory['time'].T))]
            for column in sorted_keys:
                columns.append(np.vstack((trajectory[column].T)))
            out_array = np.hstack(columns)
            out_data.append(out_array)
        return out_data

    @classmethod
    def runArray(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        self.simulation_data = []
        #create mapping of species dictionary to array indices
        species = model.listOfSpecies.keys():
        #create numpy array for timeline
        timeline = np.linspace(0, t, (t//increment+1))
        #create numpy matrix to mark other state data of species
        species_arr = np.zeros((timeline.size, len(species)))
        #copy initial values to base
        for i in range(len(species)):
            species_arr[0][i] = model.listOfSpecies[species[i]].initial_value 
        #column stack timeline and species 
        trajectoryBase = np.column_stack(timeline, species_arr)
        #create mapping of reaction dictionary to array indices
        reactions = model.listOfReactions.keys()
        propensity_functions = [None] * len(reactions)
        #pre-evaluate propensity equations from strings:
        
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        self.simulation_data = []
        curr_state = {}
        propensity = {}
        propensityFuns = {}
        for r in model.listOfReactions:
            propensityFuns[r] = eval('lambda :'+model.listOfReactions[r].propensity_function, curr_state)

        
        for traj_num in range(number_of_trajectories):
            trajectory = {}
            self.simulation_data.append(trajectory)
            trajectory['time'] = np.linspace(0,t,(t//increment+1))
            for s in model.listOfSpecies:   #Initialize Species population
                curr_state[s] = model.listOfSpecies[s].initial_value
                trajectory[s] = np.zeros(shape=(trajectory['time'].size))
            curr_state['vol'] = model.volume
            curr_time = 0	  
            entry_count = 0
            for p in model.listOfParameters:
                curr_state[p] = model.listOfParameters[p].value		

            reaction = None
            while(entry_count < trajectory['time'].size):
                prop_sum = 0
                for r in model.listOfReactions:
                    propensity[r] = (propensityFuns[r])()#eval(model.listOfReactions[r].propensity_function, curr_state)
                    prop_sum += propensity[r]
                cumil_sum = random.uniform(0,prop_sum)
                for r in model.listOfReactions:
                    cumil_sum -= propensity[r]
                    if(cumil_sum <= 0):
                        reaction = r
                        break
                if(prop_sum <= 0):
                    while(entry_count < trajectory['time'].size):
                        for s in model.listOfSpecies:
                            trajectory[s][entry_count]=(curr_state[s])
                        entry_count += 1
                    break

                tau = -1*math.log(random.random())/prop_sum
                curr_time += tau
                while(entry_count < trajectory['time'].size and curr_time >= trajectory['time'][entry_count]):
                    for s in model.listOfSpecies:
                        trajectory[s][entry_count]=(curr_state[s])
                    entry_count += 1

                for react in model.listOfReactions[reaction].reactants:
                    curr_state[str(react)] -=  model.listOfReactions[reaction].reactants[react]
                for prod in model.listOfReactions[reaction].products:
                    curr_state[str(prod)] += model.listOfReactions[reaction].products[prod]

        if show_labels:
            return self.simulation_data;
        else:
            return self.format_trajectories(self.simulation_data)
    


        
    def get_trajectories(self, outdir, debug=False, show_labels=False):
        pass
