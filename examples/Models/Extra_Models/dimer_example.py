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

# Setup the Environment
# .............................................................................
import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../')))

# MatPlotLib is used for creating custom visualizations
import matplotlib.pyplot as plt

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


# Create the Dimerization Model
# .............................................................................
def create_dimerization(parameter_values=None):
    # Initialize Model
    model = gillespy2.Model(name="Dimerization")

    # Define Variables (GillesPy2.Species)
    m = gillespy2.Species(name='monomer', initial_value=30)
    d = gillespy2.Species(name='dimer',   initial_value=0)

    # Add Variables to Model
    model.add_species([m, d])

    # Define Parameters
    k_c = gillespy2.Parameter(name='k_c', expression=0.005)
    k_d = gillespy2.Parameter(name='k_d', expression=0.08)

    # Add Parameters to Model
    model.add_parameter([k_c, k_d])

    # Define Reactions
    r_creation = gillespy2.Reaction(
        name="r_creation", reactants={'monomer': 2}, products={'dimer': 1}, rate='k_c'
    )
    r_dissociation = gillespy2.Reaction(
        name="r_dissociation", reactants={'dimer': 1}, products={'monomer': 2}, rate='k_d'
    )

    # Add Reactions to Model
    model.add_reaction([r_creation, r_dissociation])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=100, num_points=101)

    # Set Model Timespan
    model.timespan(tspan)
    return model

# Instantiate the Model
model = create_dimerization()


# Run the Simulations
# .............................................................................
# In this example, we will run the simulation 10 times, by passing the argument
# "number_of_trajectories" to the "run" method.
output = model.run(number_of_trajectories=10)


# Visualizations
# .............................................................................
# The format of the results from a model run is an array of values for
# different time points.  There will be one subarray for each trajectory.
# Here, we plot each of the 10 trajectories in the same figure.
plt.figure(figsize=(8, 5))

for trajectory in output:
    plt.plot(trajectory['time'], trajectory['monomer'], 'r')
    plt.plot(trajectory['time'], trajectory['dimer'],   'b')

# Do some additional setting up of the plot, and finally show it.

plt.legend(['M', 'D'], loc='best')
plt.xlabel("Time (s)")
plt.ylabel("Number of molecules")
plt.title("Dimerization example using GillesPy2")

plt.show()

# The end.
