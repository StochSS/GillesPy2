from gillespy2.core import GillesPySolver, Model
from gillespy2.solvers.cpp.c_simulation import CSimulation
from gillespy2.solvers.cpp.build import BuildEngine

class ODECSolver(GillesPySolver, CSimulation):
    type = "ODESimulation"

    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'timeout', 'increment', 'seed', 'debug', 'profile')

    def run(self=None, model=None, t=20, number_of_trajectories=1, timeout=0,
            increment=0.05, seed=None, debug=False, profile=False, variables={}, resume=None, **kwargs):

        pause = False