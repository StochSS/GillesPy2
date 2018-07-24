"""
An optional solver for simulating models which allows for reactions
to be effected by continuous species.

This solver utilizes RateRules in order to specify the rate
of a continuous species.  The populations of these continuous
species can in turn be incorporated into the propensity function
of one or many reactions.  This solver uses a root solving method
which sets each reaction to a uniform random number, then solves an
ODE based on the current state of the system to determine the next
fired reaction.

"""

from .gillespySolver import GillesPySolver
import random
import numpy
from scipy.integrate import ode
import math
import odespy

eval_globals = math.__dict__



class BasicHybridSolver(GillesPySolver):
    """
            Optional solver for continuous/discrete hybrid systems.

            Attributes
            ----------
            name : str
                The name by which this solver will be called.
            debug : bool
                If True, print additional debugging information to console.
            """
    name = "Basic Hybrid Solver"

    def __init__(self, debug=False):
        self.debug = debug

    @staticmethod
    def f(y, t, curr_state, reactions, rate_rules, propensities):
        """
            Evaluate the propensities for the reactions and the RHS of the RateRules.

            Attributes
            ----------
            t : float
                The time at which the function is analyzed.
            y : list
                The initial state of the system reactions and continuous species when f is called
            curr_state : dict
                Dictionary of current states for all species and parameters
            reactions : dict
                Dictionary of reactions to be evaluated
            rate_rules : dict
                Dictionary of rate_rules to be evaluated
            propensities : dict
                Dictionary of each reactions propensities to be evaluated
        """
        curr_state['t'] = t
        state_change = []
        for i, r in enumerate(reactions):
            propensities[r] = eval(reactions[r].propensity_function, eval_globals, curr_state)
            state_change.append(propensities[r])
        for i, rr in enumerate(rate_rules):
            state_change.append(eval(rate_rules[rr].expression,  eval_globals, curr_state))
        return state_change

    @staticmethod
    def reaction_fired(model):
        return lambda y,t,step_no: numpy.any( numpy.array(y[0:len(model.listOfReactions)][step_no]) > 0)
        # fired = False
        # for i, r in enumerate(model.listOfReactions):
        #     if lambda y, t, step_no: y[i][step_no] > 0:
        #         fired = True
        #fired = lambda y, t, step_no: numpy.any(numpy.array(y[0:y[-1]][step_no]) > 0)
        #return fired(y, t, step_no)

    @staticmethod
    def get_reaction_integrate(step, euler_step, curr_state, y0, model, curr_time, propensities, projected_reaction):
        """
        Helper function to perform the ODE integration of one step

        Attributes
            ----------
            step : float
                Tau distance to step during integration.
            euler_step : float
                Euler-forward predicted tau, used upon failed integration
            curr_state : dict
                Dictionary of current states for all species and parameters
            y0 : list
                The initial state of the system reactions and continuous species
            model : gillespy2.Model
                model on which the solver is run
            curr_time : float
                current time of simulation, start time of integration
            propensities : dict
                Dictionary of each reactions propensities to be evaluated
            projected_reaction : gillespy2.Reaction
                Reaction predicted by Euler-foward method, used upon failed integration
        """
        #rhs = ode(BasicHybridSolver.f)  # set function as ODE object
        #rhs.set_initial_value(y0, curr_time).set_f_params(curr_state, model.listOfReactions,
        #                                                  model.listOfRateRules, propensities)
        # rhs.set_integrator('dop853')

        #current = rhs.integrate(step+curr_time)   # current holds integration from current_time to int_time
        rhs = odespy.Vode(BasicHybridSolver.f, f_args=[curr_state, model.listOfReactions, model.listOfRateRules, propensities])
        rhs.set_initial_condition(y0)
        #time_points = numpy.arange(curr_time, euler_step, 1e-10) # make step size settable
        print("start time: ", curr_time, "step taken: ", euler_step, "number of steps: 11")
        time_points = numpy.linspace(curr_time, euler_step, 11)
        current, curr_time = rhs.solve(time_points)
        print("Current: ", current[-1, :], "\ncurr_time: ", curr_time[-1])
        return current[-1, :], curr_time[-1]
        # if rhs.successful():
        #     return current, curr_time + step
        # else:
        #     # if step is < 1e-15, take a Forward-Euler step for all species ('propensites' and RateRules)
        #     print("Projected Reaction: ", projected_reaction.name)
        #     for i, r in enumerate(model.listOfReactions):
        #         if model.listOfReactions[r] == projected_reaction:
        #             current[i] = 0
        #         else:
        #             current[i] = propensities[r] * euler_step + y0[i]
        #     # TODO The RateRule linked species should still contain the correct value in current, verify this
        #     for i, rr in enumerate(model.listOfRateRules):
        #         print("RHS FAILED: value of continuous species is: ", current[i+len(model.listOfReactions)])
        #     return current, curr_time + euler_step

    def get_reaction(self, euler_step, curr_state, y0, model, curr_time, save_time,
                     propensities, projected_reaction, debug):
        """
                Get the time to the next reaction by integrating the SSA reaction functions
            along with the RateRules.  Update population of species governed by rate rules

                Attributes
                    ----------
                    euler_step : float
                        Euler-forward predicted tau, used upon failed integration
                    curr_state : dict
                        Dictionary of current states for all species and parameters
                    y0 : list
                        The initial state of the system reactions and continuous species
                    model : gillespy2.Model
                        model on which the solver is run
                    curr_time : float
                        current time of simulation, start time of integration
                    save_time : float
                        next time point for returning data from simulation
                    propensities : dict
                        Dictionary of each reactions propensities to be evaluated
                    projected_reaction : gillespy2.Reaction
                        Reaction predicted by Euler-foward method, used upon failed integration
                    debug : bool
                        If True, print additional debugging information to console.
                """

        last_state = y0
        last_time = curr_time
        step = euler_step
        recursion_counter = 0
        time_advance_flag = False

        while True: 
            if curr_time+step > save_time:
                if debug:
                    print("Step exceeds save_time, changing step size from ", step,
                          " to ", save_time - curr_time)
                step = save_time - curr_time

            if debug:
                print("Curr Time: ", curr_time, " Save time: ", save_time,  "step: ", step)
                
            current, curr_time = self.get_reaction_integrate(step, euler_step, curr_state, y0, model,
                                                             curr_time, propensities, projected_reaction)
            print("AFTER: curr_time: ", curr_time, " current: ", current)

            occurred = []
            for i, r in enumerate(model.listOfReactions):
                if current[i] >= 0:
                    occurred.append(r)
            n_occur = len(occurred)
            if n_occur == 1:
                break
            elif n_occur > 1:
                if debug:
                    print("Multiple reactions fired in this step (n=", n_occur, ") changing step size from ", step,
                          " to ", step*0.75, "recursion_counter: ", recursion_counter)
                # reset state, and try again
                step = step * .75
                curr_time = last_time
                y0 = last_state
                recursion_counter += 1
                if recursion_counter > 20:
                    raise Exception("get_reaction() failed, too many step size reductions: halved {0} times"
                                    .format(recursion_counter))
            elif curr_time >= save_time:
                occurred.append(None)
                break
            else:
                # ODE was successful, but no reactions fired, advance time
                last_state = current
                last_time = curr_time
                if time_advance_flag:
                    if debug:
                        print("No reactions fired in this step (n=", n_occur, ") changing step size from ", step,
                              " to ", step * 1.25)
                    step = step * 1.25
                else:
                    time_advance_flag = True

        # UPDATE THE STATE of the continuous species
        for i, s in enumerate(model.listOfRateRules):
            curr_state[s] = current[i+len(model.listOfReactions)]

        if debug:
            print("Reaction Fired: ", occurred)
            print("y(t) = ", current)

        return occurred[0], current, curr_state, curr_time

    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=False,
            **kwargs):
        """
                    Runs the solver.

                    Attributes
                    ----------
                    model : Model
                        model to perform the solver on
                    t : int
                        total run time
                    number_of_trajectories : int
                        total number of trajectories to run
                    increment : float
                        time step size for along time span t
                    seed : int
                        The random seed for the simulation. Optional, defaults to None.
                    debug : bool (False)
                        Set to True to provide additional debug information about the
                        simulation.
                    show_labels : bool (True)
                        Use names of species as index of result object rather than position numbers.
                    """
        if not isinstance(self, BasicHybridSolver):
            self = BasicHybridSolver()
        if debug:
            print("t = ", t)
            print("increment = ", increment)

        random.seed(seed)

        y0 = [0] * (len(model.listOfReactions) + len(model.listOfRateRules))
        propensities = {}
        curr_state = {}
        curr_time = 0
        curr_state['vol'] = model.volume
        save_time = 0

        results = {'time': []}

        for s in model.listOfSpecies:
            # initialize populations
            curr_state[s] = model.listOfSpecies[s].initial_value
            results[s] = []

        for p in model.listOfParameters:
            curr_state[p] = model.listOfParameters[p].value

        propensity_sum = 0
        for i, r in enumerate(model.listOfReactions):
            y0[i] = (math.log(random.uniform(0, 1)))
            if debug:
                print("Setting Random number ", y0[i], " for ", model.listOfReactions[r].name)

        while save_time < t:
            while curr_time < save_time:
                projected_reaction = None
                tau_step = None
                tau_j = {}
                for i, rr in enumerate(model.listOfRateRules):
                    spec = model.listOfRateRules[rr].species.name
                    y0[i + len(model.listOfReactions)] = curr_state[spec]
                for i, r in enumerate(model.listOfReactions):
                    propensities[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
                    propensity_sum += propensities[r]
                    # Salis et al. eq (16)
                    # TODO: this needs to be optimized.  Going too big is expensive, too small is also expensive
                    if propensities[r] > 0:
                        tau_j[r] = -y0[i] / propensities[r]
                        if debug:
                            print("Propensity of ", r, " is ", propensities[r], "tau_j is ", tau_j[r])
                        if tau_step is None or tau_j[r] < tau_step:
                            tau_step = tau_j[r]
                            projected_reaction = model.listOfReactions[r]
                    else:
                        if debug:
                            print("Propensity of ", r, " is 0")
                if tau_step is None:
                    tau_step = save_time - curr_time
                if debug:
                    if projected_reaction is None:
                        print("NO projected reaction")
                    else:
                        print("Projected reaction is: ", projected_reaction.name, " at time: ", curr_time+tau_step,
                              " step size", tau_step)

                reaction, y0, curr_state, curr_time = self.get_reaction(
                   tau_step, curr_state, y0, model, curr_time, save_time, propensities, projected_reaction, debug)
                # Update curr_state with the result of the SSA reaction that fired
                if reaction is not None:
                    for i, r in enumerate(model.listOfReactions):
                        if r == reaction:
                            y0[i] = (math.log(random.uniform(0, 1)))
                            if debug:
                                print("Setting Random number ", y0[i], " for ", model.listOfReactions[r].name)
                            break
                    for reactant in model.listOfReactions[reaction].reactants:
                        curr_state[str(reactant)] -= model.listOfReactions[reaction].reactants[reactant]
                    for product in model.listOfReactions[reaction].products:
                        curr_state[str(product)] += model.listOfReactions[reaction].products[product]
            results['time'].append(save_time)
            for i, s in enumerate(model.listOfSpecies):
                results[s].append(curr_state[s])
            save_time += increment
        return results
