from .ode_solver import ODESolver


class BasicODESolver(ODESolver):
    from gillespy2.core import log
    log.warning("The name 'BasicODESolver' has been deprecated, future versions of GillesPy2 will not allow"
                " this import. Please import 'TauHybridSolver' by: from gillespy2 import ODESolver")