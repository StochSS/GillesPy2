"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from .gillespyError import SimulationError


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

    def get_increment(self, model, increment):
        """
        Set the default increment value if it was not provided

        :param model: The model on which the tspan can be found.
        :type model: gillespy.Model

        :param increment: The current value of increment.
        :type increment: int
        """
        if increment is None:
            return model.tspan[-1] - model.tspan[-2]
        if model.user_set_tspan:
            raise  SimulationError(
                """
                Failed while preparing to run the model. Both increment and timespan are set.

                To continue either remove your `timespan` definition from your Model or remove the 
                `increment` argument from this `solver.run()` call.               
                """
            )
        return increment