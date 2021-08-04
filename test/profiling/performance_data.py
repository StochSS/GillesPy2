"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

class PerformanceEntry:
    perf_time: float = 0.0
    percent: float = 0.0

    def __init__(self, t=0.0, percent=0.0):
        self.perf_time = float(t) * 1000
        self.percent = float(percent)

class PerformanceData:
    # Value in ms representing time spent overall by the executable.
    # Measures as the time from start to the time it returns.
    execution_time: float = -1.0

    # Value in ms representing time spent on the program's user-defined call stack.
    # Measured as the cumulative time spent on each sampled function of the program.
    sample_time: float = -1.0

    call_list: "dict[str, PerformanceEntry]" = {}
    worst_entry: "tuple[str, PerformanceEntry]" = ("unknown", PerformanceEntry())

    def __str__(self):
        string = f"Total execution time: {round(self.execution_time, 6)}ms ({self.sample_time}ms sampled)\n"
        for call_key, call_entry in self.call_list.items():
            string += f"[{call_key}]: {call_entry.perf_time}ms ({call_entry.percent}%)\n"

        return string
