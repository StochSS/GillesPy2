import pstats
import cProfile

from gillespy2.core import Model
from gillespy2.core import GillesPySolver
from gillespy2.solvers.cpp.c_solver import CSolver

def run_profiler(model: Model, solver: GillesPySolver, trajectories=4, timesteps=101):
    """
    Profile a Python solver running the specified Model.
    Only runtime performance is saved.
    """

    if isinstance(solver, CSolver):
        raise Exception("'solver' is actually a CSolver. This profiler cannot test C++ solvers, use 'c_profile' instead.")

    profiler = cProfile.Profile()
    profiler.enable()

    solver.run(model=model, number_of_trajectories=trajectories, t=timesteps)

    profiler.disable()

    stats = pstats.Stats(profiler).sort_stats("ncalls")
    stats.print_stats()