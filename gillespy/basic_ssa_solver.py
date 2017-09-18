import .gillespy


class BasicSSASolver(GillesPySolver):
    """ TODO
    """
    
    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False):
        """ TODO
        """
        self = BasicSSASolver()
        self.simulation_data = []
        for traj_num in range(number_of_trajectories):
            # put SSA loop here
            self.simulation_data[traj_num] = {}
            
            #to update  propensity function i   
            eval(model.listOfReactions[i].propensity_function,{},{'A':1,'B':10})
            # append
            for species_name in model.listOfSpecies:
                self.simulation_data[traj_num][species_name] = #time series data here

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if show_labels:
            return self.simulation_data
        else:
            #TODO: need to account for 'show_labels'




