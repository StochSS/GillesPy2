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
import copy
import numpy

from .timespan import TimeSpan
from .gillespyError import SimulationError, ModelError
from typing import Set, Type


class GillesPySolver:
    """ 
    Abstract class for a solver.
    """

    name = "GillesPySolver"

    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
            debug=False, profile=False, show_labels=None, **kwargs):
        """ 
        Call out and run the solver. Collect the results.

        :param model: The model on which the solver will operate.
        :type model: gillespy.Model
        
        :param t: The end time of the solver
        :type t: int
        
        :param number_of_trajectories: The number of times to sample the chemical master equation. Each
            trajectory will be returned at the end of the simulation.
        :type number_of_trajectories: int

        :param increment: The time step of the solution
        :type increment: float
    
        :param seed: The random seed for the simulation. Defaults to None. 
        :type seed: int
        
        :param debug: Set to True to provide additional debug information about the simulation.
        :type debug: bool
        
        :param show_labels: Use names of species as index of result object rather than position numbers.
        :type show_labels: bool

        :returns: Simulation trajectories.
        """
        raise SimulationError("This abstract solver class cannot be used directly.")

    def get_solver_settings(self):

        raise SimulationError("This abstract solver class cannot be used directly")

    def validate_tspan(self, increment, t):
        """
        Validate the models time span and set it if not provided.

        :param increment: The current value of increment.
        :type increment: int

        :param t: The end time of the simulation.
        :type t: int

        :raises SimulationError: if timespan and increment are both set by the user or neither are set by the user.
        """
        if self.model.tspan is None and increment is None:
            raise SimulationError(
                """
                Failed while preparing to run the model. Neither increment or timespan are set.

                To continue either add a `timespan` definition to your Model or add the 
                `increment` and `t` arguments to this `solver.run()` call.               
                """
            )

        if self.model.tspan is not None and increment is not None:
            raise  SimulationError(
                """
                Failed while preparing to run the model. Both increment and timespan are set.

                To continue either remove your `timespan` definition from your Model or remove the 
                `increment` argument from this `solver.run()` call.               
                """
            )

        if self.model.tspan is None:
            if t is None:
                tspan = TimeSpan.arange(increment)
            else:
                tspan = TimeSpan.arange(increment, t=t)
            self.model.timespan(tspan)
        elif not isinstance(self.model.tspan, TimeSpan) or type(self.model.tspan).__name__ != "TimeSpan":
            tspan = TimeSpan(self.model.tspan)
            self.model.timespan(tspan)
        else:
            self.model.tspan.validate()

    @classmethod
    def get_supported_features(cls) -> "Set[Type]":
        return set()

    @classmethod
    def validate_model(cls, sol_model, model):
        if model is not None:
            model.resolve_all_parameters()
            if model.tspan is None:
                model = copy.deepcopy(model)
                model.tspan = sol_model.tspan
            if model.get_json_hash() != sol_model.get_json_hash():
                raise SimulationError("Model must equal ODECSolver.model.")

    @classmethod
    def validate_sbml_features(cls, model):
        unsupported_features = model.get_model_features() - cls.get_supported_features()
        if unsupported_features:
            unsupported_features = [feature.__name__ for feature in unsupported_features]
            raise ModelError(f"Could not run Model, "
                             f"SBML Features not supported by {cls.name}: " +
                             ", ".join(unsupported_features))
