from .ode_solver import ODESolver


class BasicODESolver(ODESolver):
    from gillespy2.core import log
    log.warning("The name 'ODESolver' has been deprecated, future versions of GillesPy2 will not allow"
                " this import. Please import 'ODESolver'")