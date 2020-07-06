from .tau_hybrid_solver import TauHybridSolver

class BasicTauHybridSolver(TauHybridSolver):
    from gillespy2.core import log
    log.warning("The name 'BasicTauHybridSolver' has been deprecated, future versions of GillesPy2 will not allow"
                " this import. Please import 'TauHybridSolver'")
