# Module exceptions
class ModelError(Exception):
    pass
class SpeciesError(ModelError):
    pass
class ReactionError(ModelError):
    pass
class ParameterError(ModelError):
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
