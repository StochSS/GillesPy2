# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2024 GillesPy2 developers.

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

from typing import Union


class PerformanceEntry:
    perf_time: float = 0.0
    percent: float = 0.0

    def __init__(self, t=0.0, percent=0.0):
        self.perf_time = float(t) * 1000
        self.percent = float(percent)


class ExecutionData(object):
    # Value in ms representing time spent overall by the executable.
    # Measures as the time from start to the time it returns.
    execution_time: float = 0.0

    # Value in ms representing time spent on the program's user-defined call stack.
    # Measured as the cumulative time spent on each sampled function of the program.
    sample_time: float = 0.0

    def __init__(self, execution_time: "float" = 0.0, sample_time: "float" = 0.0):
        self.execution_time = execution_time
        self.sample_time = sample_time

    def __add__(self, other: "Union[int, float, ExecutionData]"):
        if isinstance(other, ExecutionData):
            return ExecutionData(
                execution_time=self.execution_time + other.execution_time,
                sample_time=self.sample_time + other.sample_time,
            )
        
        return ExecutionData(
            execution_time=self.execution_time + other,
            sample_time=self.sample_time + other,
        )

    def __truediv__(self, other: "Union[int, float, ExecutionData]"):
        if isinstance(other, ExecutionData):
            return ExecutionData(
                execution_time=self.execution_time / other.execution_time,
                sample_time=self.sample_time / other.sample_time,
            )
        
        return ExecutionData(
            execution_time=self.execution_time / other,
            sample_time=self.sample_time / other,
        )

class PerformanceData(ExecutionData):
    call_list: "dict[str, PerformanceEntry]" = None
    worst_entry: "tuple[str, PerformanceEntry]" = None

    def __init__(self, call_list=None, worst_entry=None):
        if call_list is None:
            call_list = {}
        super().__init__()

        self.call_list = call_list
        self.worst_entry = worst_entry

    def __str__(self):
        string = f"Total execution time: {round(self.execution_time, 6)}ms ({self.sample_time}ms sampled)\n"
        for call_key, call_entry in self.call_list.items():
            string += f"[{call_key}]: {call_entry.perf_time}ms ({call_entry.percent}%)\n"

        return string

class PerformanceEnsemble:
    """
    Container class for storing and operating on a sequence of profiler runs.
    Each run is made available, and statical methods can be applied on the ensemble results.
    """

    runs: "list[PerformanceData]" = None
    def __init__(self):
        self.runs = []

    def get_average_time(self) -> "ExecutionData":
        """
        Compute the average overall performance time of all profiler runs in the collection.

        :rtype: ExecutionData
        :returns: Data structure indicating the average sample time and execution time.
        """
        return self.get_total_time() / len(self.runs)

    def get_total_time(self) -> "ExecutionData":
        """
        Compute the total overall performance time of all profiler runs in the collection.

        :rtype: ExecutionData
        :returns: Data structure indicating the total sample time and execution time.
        """
        return ExecutionData(
            execution_time=sum([run.execution_time for run in self.runs]),
            sample_time=sum([run.sample_time for run in self.runs]),
        )
