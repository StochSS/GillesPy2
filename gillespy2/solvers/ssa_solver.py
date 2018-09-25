from gillespy2.solvers.gillespySolver import GillesPySolver
from gillespy2.example_models import Example
import logging

logger = logging.getLogger()

try:
    import numpy as np

    can_use_numpy = True
    from .numpy.ssa_solver import NumPySSASolver
    logging.debug("Successful Import of NumPySSASolver.")
except ModuleNotFoundError as e:
    can_use_numpy = False
    logging.warn(" Unable to use NumPy. The performance of this package can be \
                  significantly increased if you install NumPy:\nError:{0}".format(e))

try:
    import pyximport
    pyximport.install(setup_args={'include_dirs': np.get_include()})
    from .cython.cython_ssa_solver import CythonSSASolver
    can_use_cython = True
    logging.debug("Successful Import of CythonSSASolver.")
except Exception as e:
    logging.warn(" Unable to use Cython optimized SSA: {0}\nThe performance of this\
                 package can be significantly increased if you install Cython.".format(e))
    can_use_cython = False


class OptimizedSSASolver(GillesPySolver):
    """ SSA Direct Method Solver implemented primarily with Numpy. Attempts to use Cython if available.
    """
    name = "CythonSSASolver"

    def format_trajectories(simulation_data):
        out_data = []
        sorted_keys = sorted(simulation_data[0])
        sorted_keys.remove('time')
        for trajectory in simulation_data:
            columns = [np.vstack((trajectory['time'].T))]
            for column in sorted_keys:
                columns.append(np.vstack((trajectory[column].T)))
            out_array = np.hstack(columns)
            out_data.append(out_array)
        return out_data

    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=False,
            stochkit_home=None, use_cython=True):
        if use_cython and can_use_cython:
            solver = CythonSSASolver()
            return solver.run(model, t, number_of_trajectories, increment, seed, debug, show_labels)
        return super().run(model, t, number_of_trajectories, increment, seed, debug, show_labels, stochkit_home)
