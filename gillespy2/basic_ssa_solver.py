import gillespy2
import random


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

 	#reaction_list = model.get_all_reactions()
  	curr_state = {} 
	propensity = {} 
	results = {}
	#rates = {}
	#rates = dict.fromkeys(model.listOfReactions.keys())
	
        for traj_num in range(number_of_trajectories):
            # put SSA loop here

	      
            for s in model.listOfSpecies:   #Initialize Species population
                curr_state[s] = model.listOfSpecies[s].initial_value
		results[s] = []

	    curr_state['vol'] = model.volume
	    curr_time = 0
	  
            for p in model.listOfParameters:
                curr_state[p] = model.listOfParameters[p].value
	    for r in model.listOfReactions:
		print(model.listOfReactions[r].reactants)
		  
            while(curr_time < t):
		prop_sum = 0
		cumil_sum = 0
		reaction = None
		reaction_num = None
            	for r in model.listOfReactions: 
			propensity[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
			prop_sum += propensity[r] 	
		reaction_num = random.uniform(0,prop_sum)	

		for r in model.listOfReactions:
			cumil_sum += propensity[r]
			#print("Propensity: {}".format(propensity[r]))
			if(cumil_sum >= reaction_num):
				reaction = r
				#print("Reaction Num: {}".format(reaction_num))
				#print("Chosen Reaction: {}\n".format(reaction))
				break
		for reacts in model.listOfReactions[reaction].reactants:
			curr_state[reacts] -= 1
		for prods in model.listOfReactions[reaction].products:
			curr_state[prods] += 1		
			
		curr_time +=increment 
	    for s in model.listOfSpecies:
		print("{} {} ".format(s,curr_state[s])) 
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




