from .gillespyError import SimulationError


class GillesPySolver:
    name = "GillesPySolver"
    """ 
    Abstract class for a solver. This is generally called from within a
    gillespy Model through the Model.run function. Returns simulation 
    trajectories.

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
    
    :param debug: Set to True to provide additional debug information about the     
    simulation.
    :type debug: bool
    
    :param show_labels: Use names of species as index of result object rather than position numbers.
    :type show_labels: bool
    """
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
            debug=False, profile=False, show_labels=None, **kwargs):
        """ 
        Call out and run the solver. Collect the results.
        """
        raise SimulationError("This abstract solver class cannot be used directly.")

    def get_solver_settings(self):

        raise SimulationError("This abstract solver class cannot be used directly")