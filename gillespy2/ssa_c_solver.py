import gillespy2
from .gillespySolver import GillesPySolver

class SSACSolver(GillesPySolver):
    """TODO"""
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, show_labels=False,stochkit_home=None):
        #Write simulation C++ file.
        #    Write constants.
        #        Write model species names.
        #        Write initial populations.
        #        Write reaction species changes.
        #    Write propensity function.
        #Use makefile.
        #Execute simulation.
        #Parse/return results.
        pass
