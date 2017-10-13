import gillespy2
from gillespySolver import GillesPySolver
from scipy.integrate import odeint
import numpy as np
class BasicODESolver(GillesPySolver):


	@staticmethod
	def function(y0,t,species,parameters,reactions):
		curr_state = {}
   	  	state_change = {}
  		curr_state['vol'] = 1
		propensity = {}

    		for i,s in enumerate(species):  
        		curr_state[s] = y0[i]
        		state_change[s] = 0
    		for p in parameters: 
        		curr_state[p] = parameters[p].value
        
    		for r in reactions:
        		propensity[r] = eval(reactions[r].propensity_function,curr_state)  #assumption that prop is massAction
        		for react in reactions[r].reactants:
            			state_change[react] -= reactions[r].reactants[react]*propensity[r]
        		for prod in reactions[r].products:
            			state_change[prod] += reactions[r].products[prod]*propensity[r]

            	return([state_change['A'], state_change['B'], state_change['C'], state_change['D']])


	@classmethod
	def run(self, model,  t=20, number_of_trajectories=1,
		increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
		for traj_num_ in range(number_of_trajectories):
			y0 = []
			for s in model.listOfSpecies:
				y0.append(model.listOfSpecies[s].initial_value)
			time = np.arange(0,t,increment)
		return odeint(y0=y0,func=BasicODESolver.function,t=time, args=(model.listOfSpecies,model.listOfParameters,model.listOfReactions))




		
