from gillespy2.core import GillesPySolver
import random
import numpy
import math
import sys
import warnings


class BasicTauLeapingSolver(GillesPySolver):
    """
    A Basic Tau Leaping Solver for GillesPy2 models.  This solver uses an algorithm calculates
    multiple reactions in a single step over a given tau step size.  The change in propensities
    over this step are bounded by bounding the relative change in state, yielding greatly improved
    run-time performance with very little trade-off in accuracy.
    """
    name = "Basic Tau Leaping Solver"

    def __init__(self, debug=False):
        self.debug = debug
        self.epsilon = 0.03

    def __get_reactions(self, step, curr_state, curr_time, save_time, propensities, reactions, debug):
        """
        Helper Function to get reactions fired from t to t+tau.  Returns three values:
        rxn_count - dict with key=Raection channel value=number of times fired
        curr_state - dict containing all state variables for system at current time
        curr_time - float representing current time
        """


        if curr_time + step > save_time:
            if debug:
                print("Step exceeds save_time, changing step size from ", step,
                      " to ", save_time - curr_time)
            step = save_time - curr_time

        if debug:
            print("Curr Time: ", curr_time, " Save time: ", save_time, "step: ", step)

        rxn_count = {}

        for r in reactions:
            rxn_count[r] = numpy.random.poisson(propensities[r] * step)

        if debug:
            print("Reactions Fired: ", rxn_count)

        curr_time = curr_time+step

        return rxn_count, curr_state, curr_time

    def __get_tau(self, model, y0, curr_state, propensities, steps_taken, save_time, curr_time, debug, profile):
        """ Helper function to analyze best tau to take as next step """
        projected_reaction = None
        tau_step = None
        tau_j = {}

        if debug:
            print("curr_state = {", end='')
            for i, s in enumerate(model.listOfSpecies):
                print("'{0}' : {1}, ".format(s, curr_state[s]), end='')
            print("}")

        # Salis et al. eq (16)
        # TODO: this needs to be optimized.  Going too big is expensive, too small is also expensive
        propensity_sum = 0
        for i, r in enumerate(model.listOfReactions):
            propensities[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
            propensity_sum += propensities[r]
            if propensities[r] > 0:
                tau_j[r] = -y0[i] / propensities[r]
                if debug:
                    print("Propensity of ", r, " is ", propensities[r], "tau_j is ", tau_j[r])
                if tau_step is None or tau_j[r] < tau_step:
                    tau_step = max(tau_j[r], 1e-10)
                    projected_reaction = model.listOfReactions[r]
            else:
                if debug:
                    print("Propensity of ", r, " is ", propensities[r])

        if tau_step is None:
            tau_step = save_time - curr_time

        if debug:
            if projected_reaction is None:
                print("NO projected reaction")
            else:
                print("Projected reaction is: ", projected_reaction.name, " at time: ", curr_time + tau_step,
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
        for r in model.listOfReactions:
            reactant_keys = model.listOfReactions[r].reactants.keys()
            for key in reactant_keys:
                reactants.append(key)
        # initialize mean and stand_dev for reactants
        for r in reactants:
            mean[r] = 0
            stand_dev[r] = 0

        critical = False
        for r in model.listOfReactions:
            # For each reaction, determine if critical
            for reactant in model.listOfReactions[r].reactants:
                # if species pop / state change <= threshold set critical and break
                if curr_state[str(reactant)] / model.listOfReactions[r].reactants[reactant] <= n_fires:
                    critical = True
                    critical_reactions.append(r)
        if critical:
            # Cycle through critical reactions to fire fastest one, if none fire, fire soonest reaction
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
            for r in model.listOfReactions:
                for reactant in model.listOfReactions[r].reactants:
                    g_i[reactant] = 3 + (1 / (curr_state[str(reactant)] - 1)) + (
                            2 / (curr_state[str(reactant)] - 2))  # Cao, Gillespie, Petzold 27.iii
                    epsilon_i[reactant] = self.epsilon / g_i[reactant]  # Cao, Gillespie, Petzold 27
                    mean[reactant] += model.listOfReactions[r].reactants[reactant] * propensities[
                        r]  # Cao, Gillespie, Petzold 29a
                    stand_dev[reactant] += model.listOfReactions[r].reactants[reactant] ** 2 * propensities[
                        r]  # Cao, Gillespie, Petzold 29b
                for r in reactants:
                    if mean[r] > 0:
                        # Cao, Gillespie, Petzold 33
                        tau_i[r] = min((max(epsilon_i[r] * curr_state[str(r)], 1) / mean[r]),
                                       # Cao, Gillespie, Petzold 32A
                                       (max(epsilon_i[r] * curr_state[str(r)], 1) ** 2 / stand_dev[
                                           r]))  # Cao, Gillespie, Petzold 32B
                        if new_tau_step is None or tau_i[
                            r] < new_tau_step:  # set smallest tau from non-critical reactions
                            new_tau_step = tau_i[r]


        if new_tau_step is not None and new_tau_step < (
                save_time - curr_time):  # if curr+new_tau < save_time, use new_tau
            tau_step = new_tau_step
        if profile:
            steps_taken.append(tau_step)
        return new_tau_step
        # END NEW TAU SELECTION METHOD

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, profile=False,
            debug=False, show_labels=True, stochkit_home=None, **kwargs):
        """
                Function calling simulation of the model. This is typically called by the run function in GillesPy2 model
                objects and will inherit those parameters which are passed with the model as the arguments this run function.

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
            num_save_points = int(t / increment) + 1
            trajectories = numpy.empty((number_of_trajectories, num_save_points, len(model.listOfSpecies)+1))

        for trajectory in range(number_of_trajectories):
            random.seed(seed)
            y0 = [0] * (len(model.listOfReactions) + len(model.listOfRateRules))
            propensities = {}
            curr_state = {}
            curr_time = 0
            curr_state['vol'] = model.volume
            save_time = 0
            if show_labels:
                results = {'time': []}
            else:
                results = numpy.empty((number_of_trajectories, int(t / increment)+1, len(model.listOfSpecies) + 1))
            steps_taken = []
            steps_rejected = 0

            for s in model.listOfSpecies:
                # initialize populations
                curr_state[s] = model.listOfSpecies[s].initial_value
                if show_labels:
                    results[s] = []

            for p in model.listOfParameters:
                curr_state[p] = model.listOfParameters[p].value

            for i, r in enumerate(model.listOfReactions):  # set reactions to uniform random number and add to y0
                y0[i] = (math.log(random.uniform(0, 1)))
                if debug:
                    print("Setting Random number ", y0[i], " for ", model.listOfReactions[r].name)

            timestep = 0
            while save_time < t:
                while curr_time < save_time:

                    tau_step = self.__get_tau(model, y0, curr_state, propensities, steps_taken,
                                              save_time, curr_time, debug, profile)

                    prev_y0 = y0.copy()
                    prev_curr_state = curr_state.copy()
                    prev_curr_time = curr_time

                    loop_cnt = 0
                    while True:
                        loop_cnt += 1
                        if loop_cnt > 100:
                            raise Exception("Loop over __get_reactions() exceeded loop count")

                        reactions, curr_state, curr_time = self.__get_reactions(
                            tau_step, curr_state, curr_time, save_time, propensities, model.listOfReactions, debug)

                        # Update curr_state with the result of the SSA reaction that fired
                        species_modified = {}
                        for i, r in enumerate(model.listOfReactions):
                            if reactions[r] > 0:
                                for reactant in model.listOfReactions[r].reactants:
                                    species_modified[str(reactant)] = True
                                    curr_state[str(reactant)] -= model.listOfReactions[r].reactants[reactant] * reactions[r]
                                for product in model.listOfReactions[r].products:
                                    species_modified[str(product)] = True
                                    curr_state[str(product)] += model.listOfReactions[r].products[product] * reactions[r]
                        neg_state = False
                        for s in species_modified.keys():
                            if curr_state[s] < 0:
                                neg_state = True
                                if debug:
                                    print("Negative state detected: curr_state[{0}]= {1}".format(s, curr_state[s]))
                        if neg_state:
                            if debug:
                                print("\trxn={0}".format(reactions))
                            y0 = prev_y0.copy()
                            curr_state = prev_curr_state.copy()
                            curr_time = prev_curr_time
                            tau_step = tau_step / 2
                            steps_rejected += 1
                            if debug:
                                print("Resetting curr_state[{0}]= {1}".format(s, curr_state[s]))
                            if debug:
                                print("\tRejecting step, taking step of half size, tau_step={0}".format(tau_step))
                        else:
                            break  # breakout of the while True

                if show_labels:
                    results['time'].append(save_time)
                    for i, s in enumerate(model.listOfSpecies):
                        results[s].append(curr_state[s])
                else:
                    trajectories[trajectory][timestep][0] = save_time
                    for i, s in enumerate(model.listOfSpecies):
                        trajectories[trajectory][timestep][i + 1] = curr_state[s]
                save_time += increment
                timestep += 1
            if show_labels:
               trajectories.append(results)
            if profile:
                print(steps_taken)
                print("Total Steps Taken: ", len(steps_taken))
                print("Total Steps Rejected: ", steps_rejected)

        return trajectories
