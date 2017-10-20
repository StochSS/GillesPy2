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
                propensity_functions[i] = propensity_functions[i].replace(species[j], 'x[{0}]'.format(j+1))
            propensity_functions[i] = eval('lambda x:'+propensity_functions[i], parameters)
        #begin simulating each trajectory
        self.simulation_data = []
        for trajectory_num in range(number_of_trajectories):
            #copy initial state data
            trajectory = trajectory_base[trajectory_num]
            entry_count = 1
            current_time = 0
            current_state = np.copy(trajectory[0])
            propensity_sums = np.zeros(number_reactions)
            
            #calculate initial propensity sums
            while entry_count < timeline.size:
                #determine next reaction
                for i in range(number_reactions):
                    propensity_sums[i] = propensity_functions[i](current_state)
                propensity_sum = np.sum(propensity_sums)
                #if no more reactions, quit
                if propensity_sum <= 0:
                    trajectory[entry_count] = current_state
                    entry_count = timeline.size
                    break
                cumulative_sum = random.uniform(0, propensity_sum)
                current_time += -math.log(random.random()) / propensity_sum
                #determine time passed in this reaction
                while entry_count < timeline.size and timeline[entry_count] <= current_time:
                    trajectory[entry_count] = current_state
                    entry_count += 1
                for potential_reaction in range(number_reactions):
                    cumulative_sum -= propensity_sums[potential_reaction]
                    if cumulative_sum <= 0:
                            current_state[1:] += species_changes[potential_reaction]
                            #recompute propensities as needed
                            for i in range(number_reactions):
                                propensity_sums[i] = propensity_functions[i](current_state)
                            break
                
            self.simulation_data.append(trajectory)
        return self.simulation_data

        
    def get_trajectories(self, outdir, debug=False, show_labels=False):
        pass


class PrioritizedPropensity_SSASolver(GillesPySolver):
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
        
        #create numpy matrix to mark other state data of species
        species_arr = np.zeros((timeline.size, len(species)))
        #copy initial values to base
        for i in range(len(species)):
            species_arr[0][i] = model.listOfSpecies[species[i]].initial_value 
        
        trajectory_base = species_arr

        #create dictionary of all constant parameters for propensity evaluation
        parameters = {'vol' : model.volume}
        for paramName, param in model.listOfParameters.items():
            parameters[paramName] = param.value
        
        #create mapping of reaction dictionary to array indices
        reactions = list(model.listOfReactions.keys())
        number_reactions = len(reactions)
        propensity_functions = [r.propensity_function for r in model.listOfReactions.values()]

        #create an array mapping reactions to species modified
        boolean_species_changed = np.zeros((number_reactions,number_species))
        species_changes = np.zeros((number_reactions, number_species))
        #pre-evaluate propensity equations from strings:
        for i in range(number_reactions):
            #replace all references to species with array indices
            for j in range(number_species):
                if model.listOfSpecies[species[j]] in model.listOfReactions[reactions[i]].products:
                    species_changes[i][j] += model.listOfReactions[reactions[i]].products[model.listOfSpecies[species[j]]]
                    boolean_species_changed[i][j] = 1
                if model.listOfSpecies[species[j]] in model.listOfReactions[reactions[i]].reactants:
                    species_changes[i][j] -= model.listOfReactions[reactions[i]].reactants[model.listOfSpecies[species[j]]]
                    boolean_species_changed[i][j] = 1
                if species[j] in propensity_functions[i]:
                    propensity_functions[i] = propensity_functions[i].replace(species[j], 'x[{0}]'.format(j))
                    boolean_species_changed[i][j] = 1
            propensity_functions[i] = eval('lambda x:'+propensity_functions[i], parameters)
        #map reactions to species populations they change
        reaction_changes = np.zeros((number_reactions, number_reactions))
        #check all reaction changes for shared species
        for i in range(0, number_reactions):
            reaction_changes[i][i] = 1
            for j in range(i+1, number_reactions):
                if np.dot(boolean_species_changed[i], boolean_species_changed[j]) != 0:
                    reaction_changes[i][j] = 1
                    reaction_changes[j][i] = 1
        #begin simulating each trajectory
        self.simulation_data = []
        for trajectory_num in range(number_of_trajectories):
            #copy initial state data
            if trajectory_num >= number_of_trajectories - 1:
                trajectory = trajectory_base
            else:
                trajectory = np.copy(trajectory_base)

            entry_count = 1 #the initial state is 0
            current_time = 0
            current_state = np.copy(trajectory[0])

            #keep track of next reaction based on priority queue selection
            next_reaction = []
            #keep track of reactions popped off queue
            next_reaction_consumed = []
            #keep track of propensities evaluated each time step
            propensity_sums = np.zeros(number_reactions)
            #calculate initial propensity sums
            for i in range(number_reactions):
                propensity_sums[i] = propensity_functions[i](current_state)
                heapq.heappush(next_reaction, (-propensity_sums[i], i))
            while entry_count < timeline.size:
                #determine next reaction
                propensity_sum = np.sum(propensity_sums)
                #if no more reactions, quit
                if propensity_sum <= 0:
                    for entry_count in range(entry_count, timeline.size):
                        np.copyto(trajectory[entry_count], current_state)
                    break
                cumulative_sum = random.uniform(0, propensity_sum)
                current_time -= math.log(random.random()) / propensity_sum
                #determine time passed in this reaction
                while entry_count < timeline.size and timeline[entry_count] <= current_time:
                    np.copyto(trajectory[entry_count], current_state)
                    entry_count += 1
                reaction = -1
                for potential_reaction in range(number_reactions):
                    reaction_popped = heapq.heappop(next_reaction)
                    cumulative_sum += reaction_popped[0]
                    next_reaction_consumed.append(reaction_popped[1])
                    if cumulative_sum <= 0 and propensity_sums[reaction_popped[1]] > 0:
                        reaction = reaction_popped[1]
                        current_state += species_changes[reaction]
                        #recompute propensities as needed
                        for i in range(number_reactions):
                            if reaction_changes[reaction][i] != 0:
                                propensity_sums[i] = propensity_functions[i](current_state)
                        break
                #add reactions popped back to queue
                while len(next_reaction_consumed) > 0:
                    reaction = next_reaction_consumed.pop()
                    heapq.heappush(next_reaction, (-propensity_sums[reaction], reaction))
                
            self.simulation_data.append(np.column_stack((timeline,trajectory)))
        return self.simulation_data

        
    def get_trajectories(self, outdir, debug=False, show_labels=False):
        pass

