import cProfile

import pstats
from pstats import SortKey

from . import PerformanceData
from . import PerformanceEntry

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
    stats = pstats.Stats(profiler).strip_dirs().sort_stats(SortKey.TIME)

    perf_data = PerformanceData()
    perf_data.perf_time = stats.total_tt
    perf_data.call_time = stats.total_tt
    
    worst_func = stats.fcn_list[0]
    (cc, nc, tt, ct, callers) = stats.stats[worst_func]
    perf_data.worst_entry = (worst_func, PerformanceEntry(round(tt, 3), round((tt / stats.total_tt) * 100)))

    for func in stats.fcn_list:
        (cc, nc, tt, ct, callers) = stats.stats[func]

        perf_entry = PerformanceEntry(t=round(tt, 3), percent=round((tt / stats.total_tt) * 100, 3))

        parent = func[0] if func[0] != "~" else "python"
        perf_data.call_list[f"{parent}:{func[2]}"] = perf_entry
