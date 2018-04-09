import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path[:0] = ['..']
import gillespy

class parameter_changing_model(gillespy.Model):
    """
    This toy example shows how we can simply simulate the same model for 
    multiple parameter sets. Our model consists of the following reactions:

          0 -> S1
    S1 + S1 -> S2
         S2 -> 0

    So S1 being produced, dimerizing, and being degraded.
    """

    def __init__(self, parameter_values=None):

        # Initialize the model.
        gillespy.Model.__init__(self, name="simple1")
        
        # Parameters
        k1 = gillespy.Parameter(name='k1', expression=parameter_values[0])
        k2 = gillespy.Parameter(name='k2', expression=parameter_values[1])
        k3 = gillespy.Parameter(name='k3', expression=parameter_values[2])
        self.add_parameter([k1, k2, k3])
        
        # Species
        S1 = gillespy.Species(name='S1', initial_value=100)
        S2 = gillespy.Species(name='S2', initial_value=0)
        self.add_species([S1, S2])
        
        # Reactions
        rxn1 = gillespy.Reaction(
                name = 'S1 production',
                reactants = {},
                products = {S1:1},
                rate = k1 )

        rxn2 = gillespy.Reaction(
                name = 'dimer formation',
                reactants = {S1:2},
                products = {S2:1},
                rate = k2)

        rxn3 = gillespy.Reaction(
                name = 'dimer degradation',
                reactants = {S2:1},
                products = {},
                rate = k3)

        self.add_reaction([rxn1, rxn2, rxn3])
        self.timespan(np.linspace(0,100,101))



if __name__ == '__main__':

    # Here, we create the model objects. We have two different parameter 
    # sets:

    set1 = [100, 0.1, 0.1]
    set2 = [100, 0.001, 0.1]

    # For set #2, dimers (S2) form much less readily.
    
    set1_model = parameter_changing_model(parameter_values = set1)
    set2_model = parameter_changing_model(parameter_values = set2)
    
    num_trajectories = 100

    # Let's simulate for both parameter sets, and compare the results
    set1_trajectories = set1_model.run(number_of_trajectories = num_trajectories, show_labels=False)
    set2_trajectories = set2_model.run(number_of_trajectories = num_trajectories, show_labels=False)

    # Done! That was simple.









    # PLOTTING RESULTS

    # here, we will plot all trajectories with the mean overlaid
    from matplotlib import gridspec
    
    gs = gridspec.GridSpec(1,2)
    alp = 0.1 # alpha value

    # extract time values
    time = np.array(set1_trajectories[0][:,0]) 
    
    # Plot for parameter set #1
    ax0 = plt.subplot(gs[0,0])

    set1_S1 = np.array([set1_trajectories[i][:,1] for i in xrange(num_trajectories)]).T
    set1_S2 = np.array([set2_trajectories[i][:,2] for i in xrange(num_trajectories)]).T

    
    #plot individual trajectories
    ax0.plot(time, set1_S1, 'r', alpha = alp)
    ax0.plot(time, set1_S2, 'b', alpha = alp)

    #plot mean
    ax0.plot(time, set1_S1.mean(1), 'k--', label = "Mean S1")
    ax0.plot(time, set1_S2.mean(1), 'k:', label = "Mean S2")

    ax0.legend()
    ax0.set_xlabel('Time')
    ax0.set_ylabel('Species Count')
    ax0.set_title('Parameter Set 1')

    # Plot for parameter set #2
    ax1 = plt.subplot(gs[0,1])

    set2_S1 = np.array([set2_trajectories[i][:,1] for i in xrange(num_trajectories)]).T
    set2_S2 = np.array([set2_trajectories[i][:,2] for i in xrange(num_trajectories)]).T

    
    #plot individual trajectories
    ax1.plot(time, set2_S1, 'r', alpha = alp)
    ax1.plot(time, set2_S2, 'b', alpha = alp)

    #plot mean
    ax1.plot(time, set2_S1.mean(1), 'k--', label = "Mean S1")
    ax1.plot(time, set2_S2.mean(1), 'k:', label = "Mean S2")

    ax1.legend()
    ax1.set_xlabel('Time')
    ax1.set_title('Parameter Set 2')


    
    plt.tight_layout()
    plt.show()

    
    
    
    
    
    
    
    
    
    
    
    
    
