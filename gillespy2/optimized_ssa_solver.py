import gillespy2
from gillespySolver import GillesPySolver
import random
import math
import numpy as np
import heapq
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
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        #create mapping of species dictionary to array indices
        species = list(model.listOfSpecies.keys())
        number_species = len(species)

        #create numpy array for timeline
        timeline = np.linspace(0, t, (t//increment+1))
        
        #create numpy matrix to mark all state data of time and species
        trajectory_base = np.empty((number_of_trajectories,timeline.size, number_species+1))
        
        #copy time values to all trajectory row starts
        trajectory_base[:,:,0] = timeline
        #copy initial populations to base
        for i in range(number_species):
            trajectory_base[:, 0, i+1] = model.listOfSpecies[species[i]].initial_value 
        #create dictionary of all constant parameters for propensity evaluation
        parameters = {'vol' : model.volume}
        for paramName, param in model.listOfParameters.items():
            parameters[paramName] = param.value
        
        #create mapping of reaction dictionary to array indices
        reactions = list(model.listOfReactions.keys())
        number_reactions = len(reactions)
        propensity_functions = [r.propensity_function for r in model.listOfReactions.values()]
        #create an array mapping reactions to species modified
        species_changes = np.zeros((number_reactions, number_species))
        #pre-evaluate propensity equations from strings:
        for i in range(number_reactions):
            #replace all references to species with array indices
            for j in range(number_species):
                species_changes[i][j] = model.listOfReactions[reactions[i]].products.get(model.listOfSpecies[species[j]],0) - model.listOfReactions[reactions[i]].reactants.get(model.listOfSpecies[species[j]], 0)
                propensity_functions[i] = propensity_functions[i].replace(species[j], 'x[{0}]'.format(j))
            propensity_functions[i] = eval('lambda x:'+propensity_functions[i], parameters)
        #begin simulating each trajectory
        self.simulation_data = []
        for trajectory_num in range(number_of_trajectories):
            #copy initial state data
            trajectory = trajectory_base[trajectory_num]
            entry_count = 1
            current_time = 0
            current_state = np.copy(trajectory[0,1:])
            propensity_sums = np.zeros(number_reactions)
            #calculate initial propensity sums
            while entry_count < timeline.size:
                #determine next reaction
                for i in range(number_reactions):
                    propensity_sums[i] = propensity_functions[i](current_state)
                propensity_sum = np.sum(propensity_sums)
                #if no more reactions, quit
                if propensity_sum <= 0:
                    trajectory[entry_count:,1:] = current_state
                    entry_count = timeline.size
                    break
                cumulative_sum = random.uniform(0, propensity_sum)
                current_time += -math.log(random.random()) / propensity_sum
                #determine time passed in this reaction
                while entry_count < timeline.size and timeline[entry_count] <= current_time:
                    trajectory[entry_count,1:] = current_state
                    entry_count += 1
                for potential_reaction in range(number_reactions):
                    cumulative_sum -= propensity_sums[potential_reaction]
                    if cumulative_sum <= 0:
                            current_state += species_changes[potential_reaction]
                            #recompute propensities as needed
                            for i in range(number_reactions):
                                propensity_sums[i] = propensity_functions[i](current_state)
                            break
            if show_labels:
                data = {}
                data['time'] = timeline
                for i in range(number_species):
                    data[species[i]] = trajectory[:,i]
                self.simulation_data.append(data)
            else:
                self.simulation_data.append(trajectory)
        return self.simulation_data

        
    def get_trajectories(self, outdir, debug=False, show_labels=False):
        pass
