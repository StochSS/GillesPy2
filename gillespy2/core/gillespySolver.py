from .gillespyError import SimulationError


class GillesPySolver:
    name = "GillesPySolver"
    """ 
    Abstract class for a solver. This is generally called from within a
    gillespy Model through the Model.run function. Returns simulation 
    trajectories.
        
    Attributes
    ----------
    model : gillespy.Model
        The model on which the solver will operate.
    t : float
        The end time of the solver.
    number_of_trajectories : int
        The number of times to sample the chemical master equation. Each
        trajectory will be returned at the end of the simulation.
    increment : float
        The time step of the solution.
    seed : int
        The random seed for the simulation. Defaults to None.
    debug : bool (False)
        Set to True to provide additional debug information about the     
        simulation.
    show_labels : bool (True)
        Use names of species as index of result object rather than position numbers.
    """
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None,
            debug=False, profile=False, show_labels=False, **kwargs):
        """ 
        Call out and run the solver. Collect the results.
        """
        raise SimulationError("This abstract solver class cannot be used directly.")
