"""Class and methods for Basic Tau Leaping Solver"""


import random, math, sys, warnings
import signal
import numpy as np
from gillespy2.solvers.numpy import Tau
from gillespy2.core import GillesPySolver, log

from gillespy2.core.results import Results

class BasicTauLeapingSolver(GillesPySolver):
    name = 'BasicTauLeapingSolver'
    interrupted = False
    rc = 0
    """
    A Basic Tau Leaping Solver for GillesPy2 models.  This solver uses an algorithm calculates
    multiple reactions in a single step over a given tau step size.  The change in propensities
    over this step are bounded by bounding the relative change in state, yielding greatly improved
    run-time performance with very little trade-off in accuracy.
    """

    def __init__(self, debug=False, profile=False):
        name = "BasicTauLeapingSolver"
        interrupted = False
        rc = 0
        self.debug = debug
        self.profile = profile

    def get_reactions(self, step, curr_state, curr_time, save_time, propensities, reactions):
        """
        Helper Function to get reactions fired from t to t+tau.  Returns three values:
        rxn_count - dict with key=Reaction channel value=number of times fired
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
            rxn_count[rxn] = np.random.poisson(propensities[rxn] * step)

        if self.debug:
            print("Reactions Fired: ", rxn_count)

        curr_time = curr_time+step

        return rxn_count, curr_state, curr_time

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
            debug=False, profile=False, show_labels=True, display_type = None, tau_tol=0.03, **kwargs):
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
                """

        def timed_out(signum, frame):

            self.rc = 33
            self.interrupted = True

        def interval_print(signum,frame):

            import matplotlib.pyplot as plt
            from gillespy2.core.results import common_rgb_values
            from IPython.display import clear_output

            try:

                if display_type == "text":

                    # clear_output(wait=True)
                    # print("\nprogress = ", round((curr_time/timeline.size)*100,2),"%","\nt = ",str(round(curr_time,4)).ljust(8),"\n",sep="")
                    # for i in range(number_species):
                    #     print(species[i].ljust(16),":",curr_state[species[i]])

                    print(str(round(curr_time,2))[:10].ljust(10),end = "|")
                    for i in range(number_species):
                        print( str(curr_state[species[i]])[:10].ljust(10), end = "|")
                    print("")


                elif display_type == "progress":
                    clear_output(wait=True)
                    print("progress =", round((curr_time / timeline.size) * 100, 2), "%\n")

                elif display_type == "graph":

                    clear_output(wait = True)
                    plt.figure(figsize=(18,10))
                    plt.xlim(right = timeline.size)
                    for i in range(number_species):

                        line_color = common_rgb_values()[(i) % len(common_rgb_values())]

                        plt.plot(trajectory_base[0][:, 0][:entry_count].tolist() ,
                                 trajectory_base[0][:, i + 1][:entry_count].tolist(),color=line_color,label = species[i] )

                        plt.plot([entry_count- 1,curr_time],[trajectory_base[0][:, i + 1][entry_count-1] ,
                                 curr_state[species[i]]],linewidth = 3,color = line_color)

                    plt.legend(loc='upper right')
                    plt.show()

                else:
                    print("Got display_type = \"",display_type,"\". Display_type should be \"graph\", \"text\", or \"progress\"",sep="")

            except:
                print("failed to display output at curr_time =",curr_time)
                pass

        signal.signal(signal.SIGPROF, interval_print)
        signal.signal(signal.SIGALRM, timed_out)


        if not isinstance(self, BasicTauLeapingSolver):
            self = BasicTauLeapingSolver(debug=debug, profile=profile)

        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))
        if debug:
            print("t = ", t)
            print("increment = ", increment)

        species_mappings = model.sanitized_species_names()
        species = list(species_mappings.keys())
        parameter_mappings = model.sanitized_parameter_names()
        number_species = len(species)

        if seed is not None:
            if not isinstance(seed, int):
                seed = int(seed)
            if seed > 0:
                random.seed(seed)
                np.random.seed(seed)
            else:
                raise ModelError('seed must be a positive integer')

        # create numpy array for timeline
        timeline = np.linspace(0, t, int(round(t / increment + 1)))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline
        # copy initial populations to base

        for i, s in enumerate(species):
            trajectory_base[:, 0, i + 1] = model.listOfSpecies[s].initial_value
            # create dictionary of all constant parameters for propensity evaluation

        simulation_data = []

        for trajectory_num in range(number_of_trajectories):
            if self.interrupted: break
            start_state = [0] * (len(model.listOfReactions) + len(model.listOfRateRules))
            propensities = {}
            curr_state = {}
            curr_time = 0
            curr_state['vol'] = model.volume
            save_time = 0
            data = { 'time': timeline}
            steps_taken = []
            steps_rejected = 0
            entry_count = 0
            trajectory = trajectory_base[trajectory_num]


            HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, critical_threshold = Tau.initialize(model, tau_tol)

            for spec in model.listOfSpecies:
                # initialize populations
                curr_state[spec] = model.listOfSpecies[spec].initial_value

            for param in model.listOfParameters:
                curr_state[param] = model.listOfParameters[param].value

            for i, rxn in enumerate(model.listOfReactions):
                # set reactions to uniform random number and add to start_state
                start_state[i] = (math.log(random.uniform(0, 1)))
                if debug:
                    print("Setting Random number ",
                          start_state[i], " for ", model.listOfReactions[rxn].name)

            compiled_propensities = {}
            for i, r in enumerate(model.listOfReactions):
                compiled_propensities[r] = compile(model.listOfReactions[r].propensity_function, '<string>', 'eval')

            timestep = 0

            #Each save step
            while entry_count < timeline.size:
                if self.interrupted: break

                #Until save step reached
                while curr_time < save_time:
                    if self.interrupted: break
                    propensity_sum = 0

                    for i, r in enumerate(model.listOfReactions):
                        propensities[r] = eval(compiled_propensities[r], curr_state)
                        propensity_sum += propensities[r]

                    tau_args = [HOR, reactants, mu_i, sigma_i, g_i, epsilon_i, tau_tol, critical_threshold,
                            model, propensities, curr_state, curr_time, save_time]

                    tau_step = Tau.select(*tau_args)

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

                # save step reached
                for i in range(number_species):
#                     print('appending {0} to species {1}'.format(curr_state[species[i]], species[i]))
                    trajectory[entry_count][i + 1] = curr_state[species[i]]
                save_time += increment
                timestep += 1
                entry_count += 1

            # end of trajectory
            if show_labels:
                for i in range(number_species):
                    data[species[i]] = trajectory[:, i+1]
                simulation_data.append(data)
            else:
                simulation_data = trajectory_base
            if profile:
                print(steps_taken)
                print("Total Steps Taken: ", len(steps_taken))
                print("Total Steps Rejected: ", steps_rejected)

        return simulation_data, self.rc
