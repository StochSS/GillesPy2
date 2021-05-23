from gillespy2.core import GillesPySolver, Model
from gillespy2.solvers.cpp.build import BuildEngine

class ODECSolver(GillesPySolver):
    type = "ODESimulation"

    def __init__(self, model: Model = None, output_directory: str = None, delete_directory = True, resume = None, variable = True, debug = False):
        # Initialize the build engine, prepare the output directory, and compile the simulation.
        self.build_engine = BuildEngine(debug=debug, output_dir=output_directory)
        self.build_engine.prepare(model, variable=variable)
        self.executable = self.build_engine.build_simulation(self.type)

    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'timeout', 'increment', 'seed', 'debug', 'profile')

    def run(self=None, model=None, t=20, number_of_trajectories=1, timeout=0,
            increment=0.05, seed=None, debug=False, profile=False, variables={}, resume=None, **kwargs):

        pause = False