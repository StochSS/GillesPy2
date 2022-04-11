#!/usr/bin/env python3
# =============================================================================
# File:    dimer-example.py
# Summary: Reversible binding of two molecules to form a dimer
#
# -----------------------
# How to run this example
# -----------------------
#
# You can run this program from Python by starting a terminal shell on your
# computer, then changing the working directory in your shell to the
# directory where this file is located, and then running the following command:
#
#    python3 -m dimer_example.py
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
# Some molecules can self-associate to form structures known as dimers.
# "Dimerization" is a process in which two molecules of some molecular
# species (known as a "monomer" in this situation -- let's call it "M" for
# short) come together to create a new molecule (call it "D"), but do so in a
# way that is reversible, meaning the combined structure can also decay or
# dissociate back into "M". (More information about dimerization can be found
# in the online book "Molecular Biology of the Cell", 4th edition, at the
# site https://www.ncbi.nlm.nih.gov/books/NBK26830/.)  A simple model of the
# dimerization process represents it as two reactions: a reaction in which
# one molecule of "M" reacts reversibly with another molecule of "M" to form
# one new molecule (call it "D"), and another reaction in which a molecule of
# "D" breaks apart into two molecules of "M".  Each of these two reactions has
# its own rate.  In terms of biochemical reactions, it looks like this:
#
#          kc
#   2 M  <---->  D
#          kd
#
# where kc and kd represent the rate constants for creation and dissociation
# of the dimer, respectively.
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

class Dimerization(gillespy2.Model):
    def __init__(self, parameter_values=None):
        # First call the gillespy2.Model initializer.
        super().__init__(self)

        # Define parameters for the rates of creation and dissociation.
        k_c = gillespy2.Parameter(name='k_c', expression=0.005)
        k_d = gillespy2.Parameter(name='k_d', expression=0.08)
        self.add_parameter([k_c, k_d])

        # Define variables for the molecular species representing M and D.
        m = gillespy2.Species(name='monomer', initial_value=30)
        d = gillespy2.Species(name='dimer',   initial_value=0)
        self.add_species([m, d])

        # Define the reactions representing the process.  In GillesPy2,
        # the list of reactants and products for a Reaction object are each a
        # Python dictionary in which the dictionary keys are Species objects
        # and the values are stoichiometries of the species in the reaction.
        r_creation     = gillespy2.Reaction(name="r_creation", rate=k_c,
                                            reactants={m:2}, products={d:1})
        r_dissociation = gillespy2.Reaction(name="r_dissociation", rate=k_d,
                                            reactants={d:1}, products={m:2})
        self.add_reaction([r_creation, r_dissociation])

        # Set the timespan for the simulation.
        self.timespan(numpy.linspace(0, 100, 101))


# Model simulation.
# .............................................................................
# A simulation in GillesPy2 is performed by first instantiating the model to
# be simulated, and then invoking the "run" method on that object.  The
# results of the simulation are the output of "run".
#
# In this example, we will run the simulation 10 times, by passing the argument
# "number_of_trajectories" to the "run" method.

model = Dimerization()
output = model.run(number_of_trajectories=10)

# The format of the results from a model run is is an array of values for
# different time points.  There will be one subarray for each trajectory.
# Here, we plot each of the 10 trajectories in the same figure.

plt.figure(figsize=(8, 5))

for index in range(0, 10):
    trajectory = output[index]
    plt.plot(trajectory['time'], trajectory['monomer'], 'r')
    plt.plot(trajectory['time'], trajectory['dimer'],   'b')

# Do some additional setting up of the plot, and finally show it.

plt.legend(['M', 'D'], loc='best')
plt.xlabel("Time (s)")
plt.ylabel("Number of molecules")
plt.title("Dimerization example using GillesPy2")

plt.show()

# The end.
