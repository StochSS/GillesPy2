from .tau_leaping_solver import TauLeapingSolver


class BasicTauLeapingSolver(TauLeapingSolver):
    from gillespy2.core import log
    log.warning("The name 'BasicTauLeapingSolver' has been deprecated, future versions of GillesPy2 will not allow"
                " this import. Please import 'TauLeapingSolver'")