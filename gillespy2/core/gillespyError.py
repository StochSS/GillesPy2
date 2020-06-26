# Model exceptions
class ModelError(Exception):
    pass


class SpeciesError(ModelError):
    pass


class ReactionError(ModelError):
    pass


class ParameterError(ModelError):
    pass


# Solver specific errors
class SolverError(Exception):
    pass


class DirectoryError(SolverError):
    pass


class BuildError(SolverError):
    pass


class ExecutionError(SolverError):
    pass


class SimulationError(Exception):
    pass


class StochMLImportError(SimulationError):
    pass


class InvalidStochMLError(SimulationError):
    pass


class InvalidModelError(SimulationError):
    pass


class SimulationTimeoutError(SimulationError):
    pass


class EventError(ModelError):
    pass


# Results errors
class ResultsError(Exception):
    pass


class ValidationError(ResultsError):
    pass
