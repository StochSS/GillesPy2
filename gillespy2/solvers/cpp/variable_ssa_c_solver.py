from gillespy2.core import log
from . import SSACSolver


class VariableSSACSolver(SSACSolver):
    name = "VariableSSACSolver"

    def __init__(self):
        log.warning("The VariableSSACSolver has been deprecated. "
                    "Use the SSACSolver with the 'variable' parameter instead.")
