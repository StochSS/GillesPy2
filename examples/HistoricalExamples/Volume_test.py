# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:50:12 2015

@author: john
"""

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../')

import gillespy

class Simple1(gillespy.Model):
    """
    This is a simple example for mass-action degradation of species S.
    """

    def __init__(self, parameter_values=None,volume=1.0):

        # Initialize the model.
        gillespy.Model.__init__(self, name="simple1", volume=volume)
        
        # Parameters
        k1 = gillespy.Parameter(name='k1', expression=0.03)
        self.add_parameter(k1)
        
        # Species
        r1 = gillespy.Species(name='r1', initial_value=100)
        self.add_species(r1)
        r2 = gillespy.Species(name='r2', initial_value=100)
        self.add_species(r2)
        
        
        # Reactions
        rxn1 = gillespy.Reaction(
                name = 'r1d',
                reactants = {r1:2},
                products = {},
                rate = k1 )
                
        self.add_reaction(rxn1)
        
        rxn2 = gillespy.Reaction(
                name = 'r2d',
                reactants = {r2:2},
                products = {},
                propensity_function = 'k1/2 * r2*(r2-1)/vol' )
                
        self.add_reaction(rxn2)


if __name__ == '__main__':

    # Here, we create the model object.
    # We could pass new parameter values to this model here if we wished.
    simple_1 = Simple1(volume=10)

    import time
    tick = time.time()
    # The model object is simulated with the StochKit solver, and 25 
    # trajectories are returned.
    num_trajectories = 1
    simple_1trajectories = gillespy.StochKitSolver.run(simple_1, 
            number_of_trajectories = num_trajectories, show_labels=False)
    print time.time() - tick
    # PLOTTING

    # here, we will plot all trajectories with the mean overlaid
    from matplotlib import gridspec
    
    gs = gridspec.GridSpec(1,1)
    
    plt.figure()
    ax0 = plt.subplot(gs[0,0])

    # extract time values
    time = np.array(simple_1trajectories[0][:,0]) 

    # extract just the trajectories for S into a numpy array
    
    #plot mean
    ax0.plot(time,simple_1trajectories[0][:,1],'k--',label='ma')
    ax0.plot(time,simple_1trajectories[0][:,2],'g+',label='custom')
    
    ax0.legend()
    ax0.set_xlabel('Time')
    ax0.set_ylabel('Species r Count')
    
    plt.tight_layout()
    plt.show()

    
    
    
    
