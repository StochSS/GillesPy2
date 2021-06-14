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
