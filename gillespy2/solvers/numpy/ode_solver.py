"""GillesPy2 Solver for ODE solutions."""

from threading import Thread, Event
from scipy.integrate import ode
from collections import OrderedDict
import numpy as np
from gillespy2.core import GillesPySolver, log, gillespyError
from gillespy2.solvers.utilities import solverutils as nputils


class ODESolver(GillesPySolver):
    """
    This Solver produces the deterministic continuous solution via ODE.
    """
    name = "ODESolver"
    rc = 0
    stop_event = None
    result = None
    pause_event = None

    def __init__(self):
        name = "ODESolver"
        rc = 0
        stop_event = None
        pause_event = None
        result = None

    @staticmethod
    def __f(t, y, curr_state, model, c_prop):
        """
        The right hand side of the differential equation, uses scipy.integrate ode
        :param t: time as a numpy array
        :param y: species pops as a list
        :param curr_state: dictionary of eval variables
        :param model: model being simulated
        :param c_prop: precompiled reaction propensity function
        :return: integration step
        """
        curr_state[0]['t'] = t
        state_change =  OrderedDict()
        for i, species in enumerate(model.listOfSpecies):
            curr_state[0][species] = y[i]
            state_change[species] = 0
        propensity = OrderedDict()
        for r_name, reaction in model.listOfReactions.items():
            propensity[r_name] = eval(c_prop[r_name], curr_state[0])
            for react, stoich in reaction.reactants.items():
                state_change[react.name] -= propensity[r_name] * stoich
            for prod, stoich in reaction.products.items():
                state_change[prod.name] += propensity[r_name] * stoich
        state_change = list(state_change.values())
        return state_change

    @classmethod
    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'increment', 'integrator', 'integrator_options',
                'timeout')

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, show_labels=True, integrator='lsoda',
            integrator_options={}, live_output=None, live_output_options={}, timeout=None, resume=None, **kwargs):
        """

        :param model: gillespy2.model class object
        :param t: end time of simulation
        :param number_of_trajectories: Should be 1.
            This is deterministic and will always have same results
        :param increment: time step increment for plotting
        :param integrator: integrator to be used form scipy.integrate.ode. Options include 'vode', 'zvode', 'lsoda',
        'dopri5', and 'dop835'.  For more details,
        see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html

        :param integrator_options: a dictionary containing options to the scipy integrator. for a list of options,
        see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html.
            Example use: {max_step : 0, rtol : .01}

        :param timeout: If set, if simulation takes longer than timeout, will exit.
        :type timeout: int
        :param resume: Result of a previously run simulation, to be resumed
        :param live_output : str The type of output to be displayed by solver. Can be "progress", "text", or "graph".
        :param live_output_options : dictionary contains options for live_output. By default {"interval":1}.
                    "interval" specifies seconds between displaying.
                    "clear_output" specifies if display should be refreshed with each displa
        """
        if isinstance(self, type):
            self = ODESolver()
        self.stop_event = Event()
        self.pause_event = Event()

        if timeout is not None and timeout <= 0:
            timeout = None
        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))
        if number_of_trajectories > 1:
            log.warning("Generating duplicate trajectories for model with ODE Solver. "
                        "Consider running with only 1 trajectory.")

        if resume is not None:
            # start where we last left off if resuming a simulation
            lastT = resume['time'][-1]
            step = lastT - resume['time'][-2]
            timeline = np.arange(lastT, t+step, step)
        else:
            timeline = np.linspace(0, t, int(round(t / increment + 1)))

        species = list(model._listOfSpecies.keys())
        trajectory_base, tmpSpecies = nputils.numpy_trajectory_base_initialization(model, number_of_trajectories,
                                                                                   timeline, species, resume=resume)

        # curr_time and curr_state are list of len 1 so that __run receives reference
        if resume is not None:
            curr_time = [resume['time'][-1]]
        else:
            curr_time = [0]  # Current Simulation Time
        curr_state = [None]
        live_grapher = [None]

        sim_thread = Thread(target=self.___run, args=(model, curr_state, curr_time, timeline, trajectory_base,
                                                      tmpSpecies, live_grapher,), kwargs={'t': t,
                                                                                          'number_of_trajectories':
                                                                                              number_of_trajectories,
                                                                                          'increment': increment,
                                                                                          'timeout': timeout,
                                                                                          'resume': resume,
                                                                                          'integrator': integrator,
                                                                                          'integrator_options':
                                                                                              integrator_options, })
        try:
            time = 0
            sim_thread.start()
            if live_output is not None:
                import gillespy2.core.liveGraphing
                live_output_options['type'] = live_output

                gillespy2.core.liveGraphing.valid_graph_params(live_output_options)

                if resume is not None:
                    resumeTest = True  # If resuming, relay this information to live_grapher
                else:
                    resumeTest = False
                live_grapher[0] = gillespy2.core.liveGraphing.LiveDisplayer(model, timeline, number_of_trajectories,
                                                                            live_output_options, resume=resumeTest)
                display_timer = gillespy2.core.liveGraphing.RepeatTimer(live_output_options['interval'],
                                                                        live_grapher[0].display,
                                                                        args=(curr_state, curr_time, trajectory_base,live_output))
                display_timer.start()

            if timeout is not None:
                while sim_thread.is_alive():
                    sim_thread.join(.1)
                    time += .1
                    if time >= timeout:
                        break
            else:
                while sim_thread.is_alive():
                    sim_thread.join(.1)

            if live_grapher[0] is not None:
                display_timer.cancel()

            self.stop_event.set()
            while self.result is None:
                pass
        except KeyboardInterrupt:
            if live_output:
                display_timer.pause = True
                display_timer.cancel()
            self.pause_event.set()
            while self.result is None:
                pass
        if hasattr(self, 'has_raised_exception'):
            raise self.has_raised_exception
        return self.result, self.rc

    def ___run(self, model, curr_state, curr_time, timeline, trajectory_base, tmpSpecies, live_grapher, t=20,
               number_of_trajectories=1, increment=0.05, timeout=None, show_labels=True, integrator='lsoda',
               integrator_options={}, resume=None, **kwargs):
        try:
            self.__run(model, curr_state, curr_time, timeline, trajectory_base, tmpSpecies, live_grapher, t,
                       number_of_trajectories, increment, timeout, show_labels, integrator, integrator_options, resume,
                       **kwargs)
        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1

    def __run(self, model, curr_state, curr_time, timeline, trajectory_base, tmpSpecies, live_grapher, t=20,
              number_of_trajectories=1, increment=0.05, timeout=None, show_labels=True, integrator='lsoda',
              integrator_options={}, resume=None, **kwargs):

        timeStopped = 0
        if resume is not None:
            if resume[0].model != model:
                raise gillespyError.ModelError('When resuming, one must not alter the model being resumed.')
            if t < resume['time'][-1]:
                raise gillespyError.ExecutionError(
                    "'t' must be greater than previous simulations end time, or set in the run() method as the "
                    "simulations next end time")

        # compile reaction propensity functions for eval
        c_prop = OrderedDict()
        for r_name, reaction in model.listOfReactions.items():
            c_prop[r_name] = compile(reaction.ode_propensity_function, '<string>', 'eval')

        result = trajectory_base[0]
        entry_count = 0

        y0 = [0] * len(model.listOfSpecies)

        curr_state[0] = OrderedDict()

        if resume is not None:
            for i, s in enumerate(tmpSpecies):
                curr_state[0][s] = tmpSpecies[s]
                y0[i] = tmpSpecies[s]
        else:
            for i, s in enumerate(model.listOfSpecies.values()):
                curr_state[0][s.name] = s.initial_value
                y0[i] = s.initial_value

        for p_name, param in model.listOfParameters.items():
            curr_state[0][p_name] = param.value
        rhs = ode(ODESolver.__f).set_integrator(integrator, **integrator_options)
        rhs.set_initial_value(y0, curr_time[0]).set_f_params(curr_state, model, c_prop)

        while entry_count < timeline.size - 1:
            if self.stop_event.is_set():
                self.rc = 33
                break
            if self.pause_event.is_set():
                timeStopped = timeline[entry_count]
                break

            int_time = curr_time[0] + increment
            entry_count += 1
            y0 = rhs.integrate(int_time)
            curr_time[0] += increment
            for i, spec in enumerate(model.listOfSpecies):
                curr_state[0][spec] = y0[i]
                result[entry_count][i+1] = curr_state[0][spec]

        results_as_dict = {
            'time': timeline
        }
        for i, species in enumerate(model.listOfSpecies):
            results_as_dict[species] = result[:, i+1]
        results = [results_as_dict] * number_of_trajectories

        if timeStopped != 0 or resume is not None:
            results = nputils.numpy_resume(timeStopped, results, resume=resume)

        self.result = results
        return results, self.rc
