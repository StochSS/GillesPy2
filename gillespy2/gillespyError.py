# Model exceptions
class ModelError(Exception):
    pass
class SpeciesError(ModelError):
    pass
class ReactionError(ModelError):
    pass
class ParameterError(ModelError):
    pass

#Solver specific errors
class SolverError(Exception):
    pass
class DirectoryError(SolverError):
    pass
class BuildError(SolverError):
    pass
class ExecutionError(SolverError):
    pass


class SimuliationError(Exception):
    pass

# Exceptions
class StochMLImportError(Exception):
    pass

class InvalidStochMLError(Exception):
    pass
class InvalidModelError(Exception):
    pass
class SimulationError(Exception):
    pass
