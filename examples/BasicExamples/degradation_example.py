#!/usr/bin/env python3
# =============================================================================
# File:    degradation-example.py
# Summary: Simple protein decay model
#
# -----------------------
# How to run this example
# -----------------------
#
# You can run this program from Python by starting a terminal shell on your
# computer, then changing the working directory in your shell to the
# directory where this file is located, and then running the following command:
#
#    python3 -m degradation_example.py
#
# ---------------------------
# Description of this example
# ---------------------------
#
# This file contains a simple example to demonstrate the use of GillesPy2 to
# perform a stochastic simulation.  The example is biologically-motivated;
# however, in terms of biology, it is overly simplistic and does not capture
# the real-life complexity of the process being modeled -- the aim is not
# biological realism but rather to illustrate basic usage of GillesPy2.
#
# Molecules such as proteins can degrade due to various processes.  (More
# information about protein degradation can be found in the book "The Cell: A
# Molecular Approach", 2nd ed., https://www.ncbi.nlm.nih.gov/books/NBK9957/.)
# Perhaps the simplest way to model degradation is to express it as a decay
# process: every molecule of a given species of protein has some fixed
# probability of "disappearing" at any moment (which is to say, being turned
# into some other molecular species, one which we don't bother to make
# explicit in the model).  As a further simplification, we can assume the
# rate of decay is a constant function.  In terms of biochemical reactions,
# it looks like this:
#
#       k
#   P ----->
#
# where "P" represents the number of molecules of the protein species in
# question, and "k" is a constant that represents the rate of decay.
# Note that, mathematically, the model can be expressed as follows:
#
#   d P
#   ---  =  -(k * P)
#   d t
#
# where dP/dt is the rate of change of P over time.  However, for GillesPy2,
# the explicit mathematical formula does not need to be written out:
# GillesPy2 lets you expression a model directly in terms of biochemical
# reactions.  Isn't that nice?
#
# ------------------
# Author and history
# ------------------
# 2019-07-26 <mhucka@caltech.edu> Created.
# =============================================================================

import matplotlib.pyplot as plt
import numpy
import sys

try:
    import gillespy2
except:
    def error(msg): print('\033[91m' + msg + '\033[0m')
    error('┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓')
    error('┃ Could not import the gillespy2 package. Please refer to the ┃')
    error('┃ installation instructions for GillesPy2 for information on  ┃')
    error('┃ how to install the package on your computer.  Here are two  ┃')
    error('┃ location where you may find the instructions:               ┃')
    error('┃  - https://pypi.org/project/gillespy2                       ┃')
    error('┃  - https://github.com/gillesPy2/GillesPy2                   ┃')
    error('┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛')
    sys.exit()


# Model definition.
# .............................................................................
# In GillesPy2, a model to be simulated is expressed as an object having the
# parent class "Model".  Components of the model to be simulated, such as the
# reactions, molecular species, and the time span for simulation, are all
# defined within the class definition.

class ProteinDecay(gillespy2.Model):
    def __init__(self, parameter_values=None):
        # First call the gillespy2.Model initializer.
        super().__init__(self)

        # Define a parameter that represents the rate of decay.
        decayrate = gillespy2.Parameter(name='decayrate', expression=0.05)
        self.add_parameter([decayrate])

        # Define a molecular species representing the protein; "initial_value"
        # sets the number of molecules initial present of simulation.
        protein = gillespy2.Species(name='protein', initial_value=50)
        self.add_species([protein])

        # Define the reaction representing the decay process.  In GillesPy2,
        # the list of reactants and products for a Reaction object are each a
        # Python dictionary in which the dictionary keys are Species objects
        # and the values are stoichiometries of the species in the reaction.
        # (In this example, we have only one reactant, and no products.) The
        # "rate" for the Reaction object is a Parameter representing the
        # propensity of this reaction firing.)
        reaction = gillespy2.Reaction(name="decay", rate=decayrate,
                                      reactants={protein:1}, products={})
        self.add_reaction([reaction])

        # Set the timespan for the simulation.
        self.timespan(numpy.linspace(0, 100, 101))


# Model simulation.
# .............................................................................
# A simulation in GillesPy2 is performed by first instantiating the model to
# be simulated, and then invoking the "run" method on that object.  The
# results of the simulation are the output of "run".

model = ProteinDecay()
output = model.run()

# The format of the results from a model run is is an array of values for
# different time points.  The plotting code below creates a plot of the
# single trajectory found in the output.  To simplify the code slightly,
# we set the variable "trajectory" to the relevant part of the output array.

trajectory = output[0]      # There is only one species, hence one value

# Finally, we create a figure, plot the values, label the plot, and display it.

plt.figure(figsize=(8, 5))
plt.plot(trajectory['time'], trajectory['protein'], 'r', label='Protein')
plt.plot([0], [11])

plt.title("Protein decay example using GillesPy2")
plt.xlabel("Time (s)")
plt.ylabel("Number of protein molecules")
plt.legend(loc='best')

plt.show()

# The end.
