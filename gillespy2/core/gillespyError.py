# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Model exceptions
class ModelError(Exception):
    pass


class SpeciesError(ModelError):
    pass


class ParameterError(ModelError):
    pass


class ReactionError(ModelError):
    pass


class RateRuleError(ModelError):
    pass


class AssignmentRuleError(ModelError):
    pass


class EventError(ModelError):
    pass


class FunctionDefinitionError(ModelError):
    pass


class TimespanError(ModelError):
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


# Results errors
class ResultsError(Exception):
    pass


class ValidationError(ResultsError):
    pass


class SBMLError(Exception):
    pass