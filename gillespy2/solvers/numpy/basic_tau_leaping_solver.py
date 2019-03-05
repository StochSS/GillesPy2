"""Class and methods for Basic Tau Leaping Solver"""

import random
import math
import sys
import warnings
import numpy
from gillespy2.core import GillesPySolver


class BasicTauLeapingSolver(GillesPySolver):
    """
    A Basic Tau Leaping Solver for GillesPy2 models.  This solver uses an algorithm calculates
    multiple reactions in a single step over a given tau step size.  The change in propensities
    over this step are bounded by bounding the relative change in state, yielding greatly improved
    run-time performance with very little trade-off in accuracy.
    """
    name = "Basic Tau Leaping Solver"

    def __init__(self, debug=False, profile=False):
        self.debug = debug
        self.profile = profile
        self.epsilon = 0.03

    def get_reactions(self, step, curr_state, curr_time, save_time, propensities, reactions):
        """
        Helper Function to get reactions fired from t to t+tau.  Returns three values:
        rxn_count - dict with key=Raection channel value=number of times fired
        curr_state - dict containing all state variables for system at current time
        curr_time - float representing current time
        """


        if curr_time + step > save_time:
            if self.debug:
                print("Step exceeds save_time, changing step size from ", step,
                      " to ", save_time - curr_time)
            step = save_time - curr_time

        if self.debug:
            print("Curr Time: ", curr_time, " Save time: ", save_time, "step: ", step)

        rxn_count = {}

        for rxn in reactions:
            rxn_count[rxn] = numpy.random.poisson(propensities[rxn] * step)

        if self.debug:
            print("Reactions Fired: ", rxn_count)

        curr_time = curr_time+step

        return rxn_count, curr_state, curr_time

    def get_tau(self, model, start_state, curr_state,
                propensities, steps_taken, save_time, curr_time):
        """ Helper function to analyze best tau to take as next step """
        #   pylint: disable=W0123
        projected_reaction = None
        tau_step = None
        tau_j = {}

        if self.debug:
            print("curr_state = {", end='')
            for i, spec in enumerate(model.listOfSpecies):
                print("'{0}' : {1}, ".format(spec, curr_state[spec]), end='')
            print("}")

        # Salis et al. eq (16)
        propensity_sum = 0
        for i, rxn in enumerate(model.listOfReactions):
            propensities[rxn] = eval(model.listOfReactions[rxn].propensity_function, curr_state)
            propensity_sum += propensities[rxn]
            if propensities[rxn] > 0:
                tau_j[rxn] = -start_state[i] / propensities[rxn]
                if self.debug:
                    print("Propensity of ", rxn, " is ", propensities[rxn], "tau_j is ", tau_j[rxn])
                if tau_step is None or tau_j[rxn] < tau_step:
                    tau_step = max(tau_j[rxn], 1e-10)
                    projected_reaction = model.listOfReactions[rxn]
            else:
                if self.debug:
                    print("Propensity of ", rxn, " is ", propensities[rxn])

        if tau_step is None:
            tau_step = save_time - curr_time

        if self.debug:
            if projected_reaction is None:
                print("NO projected reaction")
            else:
                print("Projected reaction is: ",
                      projected_reaction.name, " at time: ", curr_time + tau_step,
                      " step size: ", tau_step)

        # BEGIN NEW TAU SELECTION METHOD
        g_i = {}  # used for relative error calculation
        epsilon_i = {}  # relative error allowance of species
        tau_i = {}  # estimated tau based on depletion of species
        reactants = []  # a list of all species in the model which act as reactants
        mean = {}  # mu_i for each species
        stand_dev = {}  # sigma_i squared for each species
        critical_reactions = []
        new_tau_step = None
        n_fires = 2  # if a reaction would deplete a resource in n_fires, it is considered critical

        # Create list of all reactants
        for rxn in model.listOfReactions:
            reactant_keys = model.listOfReactions[rxn].reactants.keys()
            for key in reactant_keys:
                reactants.append(key)
        # initialize mean and stand_dev for reactants
        for reactant in reactants:
            mean[reactant] = 0
            stand_dev[reactant] = 0

        critical = False
        for rxn in model.listOfReactions:
            # For each reaction, determine if critical
            for reactant in model.listOfReactions[rxn].reactants:
                # if species pop / state change <= threshold set critical and break
                if curr_state[str(reactant)] / model.listOfReactions[
                        rxn].reactants[reactant] <= n_fires:
                    critical = True
                    critical_reactions.append(rxn)
        if critical:
            # Cycle through critical reactions to fire fastest one,
            # if none fire, fire soonest reaction
            for reaction in critical_reactions:
                if propensities[reaction] > 0:
                    if new_tau_step is None:
                        new_tau_step = tau_j[reaction]
                    else:
                        if tau_j[reaction] < new_tau_step:
                            new_tau_step = tau_j[reaction]
            if new_tau_step is None:
                new_tau_step = tau_step
        else:
            for rxn in model.listOfReactions:
                for reactant in model.listOfReactions[rxn].reactants:
                    g_i[reactant] = 3 + (1 / (curr_state[str(reactant)] - 1)) + (
                        2 / (curr_state[str(reactant)] - 2))  # Cao, Gillespie, Petzold 27.iii
                    epsilon_i[reactant] = self.epsilon / g_i[reactant]  # Cao, Gillespie, Petzold 27
                    mean[reactant] += model.listOfReactions[rxn].reactants[reactant] * propensities[
                        rxn]  # Cao, Gillespie, Petzold 29a
                    stand_dev[reactant] += model.listOfReactions[rxn].reactants[
                        reactant] ** 2 * propensities[rxn]  # Cao, Gillespie, Petzold 29b
                for reactant in reactants:
                    if mean[reactant] > 0:
                        # Cao, Gillespie, Petzold 33
                        tau_i[reactant] = min((max(
                            epsilon_i[reactant] * curr_state[str(reactant)], 1) / mean[reactant]),
                                              # Cao, Gillespie, Petzold 32A
                                              (max(epsilon_i[reactant] * curr_state[
                                                  str(reactant)], 1) ** 2 / stand_dev[
                                                      reactant]))  # Cao, Gillespie, Petzold 32B
                        if new_tau_step is None or tau_i[
                                reactant] < new_tau_step:
                            # set smallest tau from non-critical reactions
                            new_tau_step = tau_i[reactant]


        if new_tau_step is not None and new_tau_step < (
                save_time - curr_time):  # if curr+new_tau < save_time, use new_tau
            tau_step = new_tau_step
        if self.profile:
            steps_taken.append(tau_step)
        return new_tau_step
        # END NEW TAU SELECTION METHOD

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
            debug=False, profile=False, show_labels=True, stochkit_home=None, **kwargs):
        """
        Function calling simulation of the model.
        This is typically called by the run function in GillesPy2 model objects
        and will inherit those parameters which are passed with the model
        as the arguments this run function.

                Attributes
                ----------

                model : GillesPy2.Model
                    GillesPy2 model object to simulate
                t : int
                    Simulation run time
                number_of_trajectories : int
                    The number of times to sample the chemical master equation. Each
                    trajectory will be returned at the end of the simulation.
                    Optional, defaults to 1.
                increment : float
                    Save point increment for recording data
                seed : int
                    The random seed for the simulation. Optional, defaults to None.
                debug : bool (False)
                    Set to True to provide additional debug information about the
                    simulation.
                profile : bool (Fasle)
                    Set to True to provide information about step size (tau) taken at each step.
                show_labels : bool (True)
                    Use names of species as index of result object rather than position numbers.
                stochkit_home : str
                    Path to stochkit. This is set automatically upon installation, but
                    may be overwritten if desired.
                """
        if not sys.warnoptions:
            warnings.simplefilter("ignore")
        if not isinstance(self, BasicTauLeapingSolver):
            self = BasicTauLeapingSolver()
        if debug:
            print("t = ", t)
            print("increment = ", increment)

        if show_labels:
            trajectories = []
        else:
            num_save_points = int(t / increment)
            trajectories = numpy.empty((number_of_trajectories,
                                        num_save_points, len(model.listOfSpecies)+1))

        for trajectory in range(number_of_trajectories):
            random.seed(seed)
            start_state = [0] * (len(model.listOfReactions) + len(model.listOfRateRules))
            propensities = {}
            curr_state = {}
            curr_time = 0
            curr_state['vol'] = model.volume
            save_time = 0
            if show_labels:
                results = {'time': []}
            else:
                results = numpy.empty((number_of_trajectories, int(t / increment)+1,
                                       len(model.listOfSpecies) + 1))
            steps_taken = []
            steps_rejected = 0

            for spec in model.listOfSpecies:
                # initialize populations
                curr_state[spec] = model.listOfSpecies[spec].initial_value
                if show_labels:
                    results[spec] = []

            for param in model.listOfParameters:
                curr_state[param] = model.listOfParameters[param].value

            for i, rxn in enumerate(model.listOfReactions):
                # set reactions to uniform random number and add to start_state
                start_state[i] = (math.log(random.uniform(0, 1)))
                if debug:
                    print("Setting Random number ",
                          start_state[i], " for ", model.listOfReactions[rxn].name)

            timestep = 0
            while save_time < t:
                while curr_time < save_time:

                    tau_step = self.get_tau(
                        model, start_state, curr_state, propensities, steps_taken,
                        save_time, curr_time)

                    prev_start_state = start_state.copy()
                    prev_curr_state = curr_state.copy()
                    prev_curr_time = curr_time

                    loop_cnt = 0
                    while True:
                        loop_cnt += 1
                        if loop_cnt > 100:
                            raise Exception("Loop over get_reactions() exceeded loop count")

                        reactions, curr_state, curr_time = self.get_reactions(
                            tau_step, curr_state, curr_time, save_time,
                            propensities, model.listOfReactions)

                        # Update curr_state with the result of the SSA reaction that fired
                        species_modified = {}
                        for i, rxn in enumerate(model.listOfReactions):
                            if reactions[rxn] > 0:
                                for reactant in model.listOfReactions[rxn].reactants:
                                    species_modified[str(reactant)] = True
                                    curr_state[str(reactant)] -= model.listOfReactions[
                                        rxn].reactants[reactant] * reactions[rxn]
                                for product in model.listOfReactions[rxn].products:
                                    species_modified[str(product)] = True
                                    curr_state[str(product)] += model.listOfReactions[
                                        rxn].products[product] * reactions[rxn]
                        neg_state = False
                        for spec in species_modified:
                            if curr_state[spec] < 0:
                                neg_state = True
                                if debug:
                                    print("Negative state detected: curr_state[{0}]= {1}".format(
                                        spec, curr_state[spec]))
                        if neg_state:
                            if debug:
                                print("\trxn={0}".format(reactions))
                            start_state = prev_start_state.copy()
                            curr_state = prev_curr_state.copy()
                            curr_time = prev_curr_time
                            tau_step = tau_step / 2
                            steps_rejected += 1
                            if debug:
                                print("Resetting curr_state[{0}]= {1}".format(
                                    spec, curr_state[spec]))
                            if debug:
                                print(
                                    "\tRejecting step, taking step of half size, ",
                                    "tau_step={0}".format(tau_step))
                        else:
                            break  # breakout of the while True

                if show_labels:
                    results['time'].append(save_time)
                    for i, spec in enumerate(model.listOfSpecies):
                        results[spec].append(curr_state[spec])
                else:
                    trajectories[trajectory][timestep][0] = save_time
                    for i, spec in enumerate(model.listOfSpecies):
                        trajectories[trajectory][timestep][i + 1] = curr_state[spec]
                save_time += increment
                timestep += 1
            if show_labels:
                trajectories.append(results)
            if profile:
                print(steps_taken)
                print("Total Steps Taken: ", len(steps_taken))
                print("Total Steps Rejected: ", steps_rejected)

        return trajectories
