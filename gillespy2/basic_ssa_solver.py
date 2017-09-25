import gillespy2


class BasicSSASolver(gillespy2.GillesPySolver):
    """ TODO
    """
    
    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        """ TODO
        """
        self = BasicSSASolver()
        self.simulation_data = []

 	reaction_list = model.get_all_reactions()
  	curr_state = {} 
	propensity = {} 
	rates = {}
	rates = dict.fromkeys(model.listOfReactions.keys())
	print(rates)
        for traj_num in range(number_of_trajectories):
            # put SSA loop here

	       
            for s in model.listOfSpecies:   #Initialize Species population
                curr_state[s] = model.listOfSpecies[s].initial_value
	       #print(model.listOfSpecies[s].initial_value)
	    curr_state['vol'] = model.volume
            for p in model.listOfParameters:
		#rates[p] = model.listOfParameters[p].value
                curr_state[p] = model.listOfParameters[p].value        
		#curr_state[s] = model.listOfParameters[p].value  
                		
            for r in model.listOfReactions:
		#propensity[r] =
		#print(model.listOfReactions[r]) 
		propensity[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
	    return propensity	
		#print(propensity[r])
		#propensity.append(eval(model.listOfReactions[r].propensity_function(), curr_state))
           
		 #self.simulation_data[traj_num] = {}
            
            #to update  propensity function i   
          #  eval(model.listOfReactions[i].propensity_function,{},{'A':1,'B':10})
            # append
         #   for species_name in model.listOfSpecies:
 #               self.simulation_data[traj_num][species_name] = #time series data here

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if show_labels:
            return self.simulation_data
       # else:
            #TODO: need to account for 'show_labels'




