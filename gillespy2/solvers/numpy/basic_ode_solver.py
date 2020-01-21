"""GillesPy2 Solver for ODE solutions."""

import signal
from scipy.integrate import ode
from scipy.integrate import odeint
from collections import OrderedDict
import numpy as np
from gillespy2.core import GillesPySolver, log


class BasicODESolver(GillesPySolver):
    """
    Solver to produce a deterministic continuous solution using an ordinary
    differential equation (ODE) solver.
    """
    name = "BasicODESolver"

    def __init__(self):
        name = "BasicODESolver"
        interrupted = False
        rc = 0

    @staticmethod
    def __f(t, y, curr_state, model, ode_species, prop_functions):
        """
        The right hand side of the differential equation. This uses
        odeint from scipy.integrate.

        :param t: time as a numpy array
        :param y: species populations as a list
        :param current_state: dictionary of eval variables
        :param model: model being simulated
        :param ode_species: list of non-boundary, non-constant species
        :param prop_functions: precompiled reaction propensity function

        :return: integration step
        """
        curr_state['t'] = t
        state_change =  OrderedDict()
        for i, species in enumerate(ode_species):
            curr_state[species.name] = y[i]
            state_change[species.name] = 0
        propensity = OrderedDict()
        for reaction in model.listOfReactions.values():
            propensity[reaction.name] = eval(prop_functions[reaction.name], curr_state)
            for reactant, stoich in reaction.reactants.items():
                if reactant.name in state_change:
                    state_change[reactant.name] -= propensity[reaction.name] * stoich
            for product, stoich in reaction.products.items():
                if product.name in state_change:
                    state_change[product.name] += propensity[reaction.name] * stoich
        state_change = list(state_change.values())
        return state_change

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05,
            show_labels=True, integrator='lsoda', integrator_options={}, **kwargs):
        """
        Run the ODE solver, and collect the results.

        :param model: gillespy2.model class object
        :param t: end time of simulation
        :param number_of_trajectories: Should be 1.
            This is deterministic and will always have same results
        :param increment: time step increment for plotting
        :param show_labels: If true, simulation returns a list of trajectories,
            where each list entry is a dictionary containing key value pairs of
            species : trajectory.  If false, returns a numpy array with
            shape [traj_no, time, species].
        :param integrator: integrator to be used form scipy.integrate.ode.
            Options include 'vode', 'zvode', 'lsoda', 'dopri5', and 'dop835'.
            For more details, see
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
        :param integrator_options: a dictionary containing options to the
            scipy integrator. For a list of options,
            see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html.
            Example use: {max_step : 0, rtol : .01}
        :param kwargs:

        :return:
        """

        # Preliminary setup and sanity checks.
        if not isinstance(self, BasicODESolver):
            self = BasicODESolver()
        if number_of_trajectories > 1:
            log.warning("Generating duplicate trajectories for model with"
                        "ODE solver. Consider running with only 1 trajectory.")
        for key in kwargs:  # Only executes if len(kwargs) > 0.
            log.warning('Unsupported keyword argument to {0}: {1}'.format(self.name, key))

        # The final output needs to contain all the species, but we don't
        # construct ODEs for species marked as boundary or constant.
        # Thus, during computation of ODE values, we use a subset of species.
        # These are the only ones that will appear on the LHS of the ODEs.
        ode_species = [s for s in model.listOfSpecies.values()
                       if not (s.is_boundary or s.is_constant)]

        # Create numpy array for the timeline.
        timeline = np.linspace(0, t, round(t / increment + 1))

        # Create a numpy matrix to mark all state data of time and species and
        # initialize the first column with the timeline values.
        trajectory = np.empty((timeline.size, len(ode_species) + 1))
        trajectory[:, 0] = timeline

        # Initialize the trajectory base with species populations.
        # Start numbering by 1, because index 0 is the time point value.
        # Also record the column index of each species for use at the end.
        # (We do that by just adding a property to the object for our purposes.)
        for i, species in enumerate(ode_species, 1):
            trajectory[0, i] = species.initial_value
            species.trajectory_index = i

        # y0 is a vector of initial values.  curr_state contains the current
        # state values and is handed as the context in each function
        # evaluation call in the integrator function __f().  The values of
        # curr_state include the species in the ODEs, the parameters, and
        # constant or boundary species that do not appear in the ODEs but can
        # appear in the RHS of the ODEs.
        y0 = []
        curr_time = 0
        curr_state = OrderedDict()
        for species in ode_species:
            curr_state[species.name] = species.initial_value
            y0.append(species.initial_value)
        for species in set(model.listOfSpecies.values()) - set(ode_species):
            curr_state[species.name] = species.initial_value
        for param in model.listOfParameters.values():
            curr_state[param.name] = param.value

        # Compile reaction propensity functions to interpretable code.
        prop_functions = OrderedDict()
        for r in model.listOfReactions.values():
            prop_functions[r.name] = compile(r.ode_propensity_function, '<string>', 'eval')

        # Pull all the pieces together to set up the ODE solver, but only if
        # we have something to integrate.  If not, skip it.
        if len(ode_species) > 0:
            rhs = ode(BasicODESolver.__f)
            rhs.set_integrator(integrator, **integrator_options)
            rhs.set_f_params(curr_state, model, ode_species, prop_functions)
            rhs.set_initial_value(y0, curr_time)

            # Now run the integrator to compute values for every trajectory point.
            for entry_count in range(1, timeline.size):
                int_time = curr_time + increment
                y0 = rhs.integrate(int_time)
                curr_time += increment
                for i, species in enumerate(ode_species):
                    curr_state[species.name] = y0[i]
                    trajectory[entry_count][i+1] = curr_state[species.name]

        # Produce results in the desired form.  Note we're careful to output
        # all species, not only the ones that are non-constant, non-boundary.
        len_timeline = len(timeline)
        if show_labels:
            results_dict = { 'time': timeline }
            # We iterate over the original species list to get correct order.
            for i, species in enumerate(model.listOfSpecies.values()):
                if species in ode_species:
                    results_dict[species.name] = trajectory[:, species.trajectory_index]
                else:
                    results_dict[species.name] = [curr_state[species.name]] * len_timeline
            return [results_dict] * number_of_trajectories
        else:
            results = np.empty((timeline.size, len(model.listOfSpecies) + 1))
            results[:, 0] = timeline
            for i, species in enumerate(model.listOfSpecies.values(), 1):
                if species in ode_species:
                    results[:, i] = trajectory[:, species.trajectory_index]
                else:
                    results[:, i] = [curr_state[species.name]] * len_timeline

            return np.stack([results] * number_of_trajectories, axis=0)
