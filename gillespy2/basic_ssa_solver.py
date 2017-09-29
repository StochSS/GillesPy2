import gillespy2
import random
import math

class BasicSSASolver(gillespy2.GillesPySolver):
    """ TODO
    """
    
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        """ TODO
        """
       # self = BasicSSASolver()
        self.simulation_data = []
 
  	curr_state = {} 
	propensity = {} 
	results = {}	
	
        for traj_num in range(number_of_trajectories):
            # put SSA loop here
	    
	      
            for s in model.listOfSpecies:   #Initialize Species population
                curr_state[s] = model.listOfSpecies[s].initial_value	
		results[s]=[]

	    curr_state['vol'] = model.volume
	    results['time'] = []
	    curr_time = 0
	    save_time = 0	  

            for p in model.listOfParameters:
                curr_state[p] = model.listOfParameters[p].value		
		
	     
            while(curr_time < t):
		prop_sum = 0
		cumil_sum = 0
		reaction = None
		reaction_num = None
            	for r in model.listOfReactions: 
			propensity[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
			prop_sum += propensity[r]
			#print('------------------------------------------propensity-------------------------') 
			#print('time',curr_time,'r',r,propensity[r], model.listOfReactions[r].propensity_function)	
		reaction_num = random.uniform(0,prop_sum)
		#print('Reaction Rand: ',reaction_num)	
		#print('Propensity Sum: ',prop_sum)
		#for r in model.listOfReactions:
		#	print(r,propensity[r])
		for r in model.listOfReactions:
			cumil_sum += propensity[r]
	
			if(cumil_sum >= reaction_num):
				reaction = r
				#print('Cumilative Sum:  ',cumil_sum)	
				break

		#print('Reaction: ',reaction)
		#print('-'*80)
	
		if(prop_sum <= 0):
			while(save_time <= t):
				results['time'].append(save_time)
				for s in model.listOfSpecies:
					results[s].append(curr_state[s])
	                        save_time += increment
	
			return results
		

		tau = -1*math.log(random.random())/prop_sum
		curr_time += tau
		while(curr_time > save_time and curr_time <= t):
			results['time'].append(save_time)
			for s in model.listOfSpecies:
				results[s].append(curr_state[s])
                        save_time += increment	

		for react in model.listOfReactions[reaction].reactants:
			curr_state[react] -=  model.listOfReactions[reaction].reactants[react]
		for prod in model.listOfReactions[reaction].products:
			curr_state[prod] += model.listOfReactions[reaction].products[prod] 	


	    return results
		
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




