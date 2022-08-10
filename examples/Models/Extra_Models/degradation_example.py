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


# Create the Protein Decay Model
# .............................................................................
def create_protein_decay(parameter_values=None):
    # Instantiate Model
    model = gillespy2.Model(name="Protein Decay")

    # Define Variables (GillesPy2.Species)
    protein = gillespy2.Species(name='protein', initial_value=50)

    # Add Variables to Model
    model.add_species(protein)

    # Define Parameters
    decayrate = gillespy2.Parameter(name='decayrate', expression=0.05)

    # Add Parameters to Model
    model.add_parameter(decayrate)

    # Define Reactions
    reaction = gillespy2.Reaction(
        name="decay", rate='decayrate', reactants={'protein': 1}, products={}
    )

    # Add Reactions to Model
    model.add_reaction(reaction)

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=100, num_points=101)

    # Set Model Timespan
    model.timespan(tspan)
    return model

# Instantiate the Model
model = create_protein_decay()


# Run the Simulations
# .............................................................................
trajectory = model.run()


# Visualizations
# .............................................................................
plt.figure(figsize=(8, 5))
plt.plot(trajectory['time'], trajectory['protein'], 'r', label='Protein')
plt.plot([0], [11])

plt.title("Protein decay example using GillesPy2")
plt.xlabel("Time (s)")
plt.ylabel("Number of protein molecules")
plt.legend(loc='best')

plt.show()

# The end.
