# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import numpy as np

from gillespy2.solvers.cpp.c_decoder import IterativeSimDecoder
from gillespy2.solvers.utilities import solverutils as cutils
from gillespy2.core import GillesPySolver, Model
from gillespy2.core.gillespyError import *
from gillespy2.core import Results 

from .c_solver import CSolver, SimulationReturnCode

class TauLeapingCSolver(GillesPySolver, CSolver):

    """
    A Tau Leaping solver for GillesPy2 models.  This solver uses an algorithm that calculates
    multiple reactions in a single step over a given tau step size.  The change in propensities
    over this step are bounded by relative change in state, yielding greatly improved
    run-time performance with very little trade-off in accuracy.
    """

    name = "TauLeapingCSolver"
    target = "tau_leap"

    @classmethod
    def get_solver_settings(cls):
        """
        Returns a list of arguments supported by tau_leaping_c_solver.run.
        :returns: Tuple of strings, denoting all keyword argument for this solvers run() method.
        :rtype: tuple
        """
        return ('model', 't', 'number_of_trajectories', 'timeout', 'increment', 'seed', 'debug', 'profile')

    def run(self=None, model: Model = None, t: int = None, number_of_trajectories: int = 1, timeout: int = 0,
            increment: int = None, seed: int = None, debug: bool = False, profile: bool = False, variables={},
            resume=None, live_output: str = None, live_output_options: dict = {}, tau_tol=0.03, **kwargs):

        """
        :param model: The model on which the solver will operate. (Deprecated)
        :type model: gillespy2.Model

        :param t: End time of simulation.
        :type t: int

        :param number_of_trajectories: Number of trajectories to simulate. By default number_of_trajectories = 1.
        :type number_of_trajectories: int

        :param timeout: If set, if simulation takes longer than timeout, will exit.
        :type timeout: int

        :param increment: Time step increment for plotting.
        :type increment: float

        :param seed: The random seed for the simulation. Optional, defaults to None.
        :type seed: int

        :param variables: Dictionary of species and their data that will override existing species data.
        :type variables: dict

        :param resume: Result of a previously run simulation, to be resumed.
        :type resume: gillespy2.Results

        :param live_output: The type of output to be displayed by solver. Can be "progress", "text", or "graph".
        :type live_output: str

        :param live_output_options: dictionary contains options for live_output. By default {"interval":1}.
            "interval" specifies seconds between displaying.
            "clear_output" specifies if display should be refreshed with each display.
        :type live_output_options:  dict

        :param tau_tol: Tolerance level for Tau leaping algorithm.  Larger tolerance values will
        result in larger tau steps. Default value is 0.03.
        :type tau_tol: float

        :returns: A result object containing the results of the simulation
        :rtype: gillespy2.Results
        """

        if self is None:
            # Post deprecation block
            # raise SimulationError("TauLeapingCSolver must be instantiated to run the simulation")
            # Pre deprecation block
            log.warning(
                """
                `gillespy2.Model.run(solver=TauLeapingCSolver)` is deprecated.

                You should use `gillespy2.Model.run(solver=TauLeapingCSolver(model=gillespy2.Model))
                Future releases of GillesPy2 may not support this feature.
                """
            )
            self = TauLeapingCSolver(model, resume=resume)

        if model is not None:
            log.warning('model = gillespy2.model is deprecated. Future releases '
                        'of GillesPy2 may not support this feature.')
        if self.model is None:
            if model is None:
                raise SimulationError("A model is required to run the simulation.")
            self._set_model(model=model)

        self.model.compile_prep()
        self.validate_model(self.model, model)
        self.validate_sbml_features(model=self.model)

        self.validate_tspan(increment=increment, t=t)
        if increment is None:
            increment = self.model.tspan[-1] - self.model.tspan[-2]
        if t is None:
            t = self.model.tspan[-1]

        # Validate parameters prior to running the model.
        self._validate_type(variables, dict, "'variables' argument must be a dictionary.")
        self._validate_variables_in_set(variables, self.species + self.parameters)
        self._validate_resume(t, resume)
        self._validate_kwargs(**kwargs)

        if resume is not None:
            t = abs(t - int(resume["time"][-1]))

        number_timesteps = int(round(t / increment + 1))

        args = {
            "trajectories": number_of_trajectories,
            "timesteps": number_timesteps,
            "tau_tol": tau_tol,
            "end": t,
            "interval": str(number_timesteps),
        }

        if self.variable:
            populations = cutils.update_species_init_values(self.model.listOfSpecies, self.species, variables, resume)
            parameter_values = cutils.change_param_values(self.model.listOfParameters, self.parameters, self.model.volume, variables)

            args.update({
                "init_pop": populations,
                "parameters": parameter_values
            })

        seed = self._validate_seed(seed)
        if seed is not None:
            args.update({
                "seed": seed
            })

        if live_output is not None:
            live_output_options['type'] = live_output
            display_args = {
                "model": self.model, "number_of_trajectories": number_of_trajectories, "timeline": np.linspace(0, t, number_timesteps),
                "live_output_options": live_output_options, "resume": bool(resume)
            }
        else:
            display_args = None

        args = self._make_args(args)
        decoder = IterativeSimDecoder.create_default(number_of_trajectories, number_timesteps, len(self.model.listOfSpecies))

        sim_exec = self._build(self.model, self.target, self.variable, False)
        sim_status = self._run(sim_exec, args, decoder, timeout, display_args)

        if sim_status == SimulationReturnCode.FAILED:
            raise ExecutionError("Error encountered while running simulation C++ file:\n"
                f"Return code: {int(sim_status)}.\n")

        trajectories, time_stopped = decoder.get_output()

        simulation_data = self._format_output(trajectories)
        if sim_status == SimulationReturnCode.PAUSED:
            simulation_data = self._make_resume_data(time_stopped, simulation_data, t)
        if resume is not None:
            simulation_data = self._update_resume_data(resume, simulation_data, time_stopped)
        self.result = simulation_data
        self.rc = int(sim_status)

        return Results.build_from_solver_results(self, live_output_options)
