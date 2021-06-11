class PerformanceEntry:
    perf_time: float = 0.0
    percent: float = 0.0

    def __init__(self, t=0.0, percent=0.0):
        self.perf_time = float(t) * 1000
        self.percent = float(percent)

class PerformanceData:
    perf_time: float = -1.0
    call_time: float = -1.0
    call_list: "dict[str, PerformanceEntry]" = {}
    worst_entry: "tuple[str, PerformanceEntry]" = ("unknown", PerformanceEntry())

    def __str__(self):
        string = f"Total time: {round(self.perf_time, 6)}ms ({self.call_time}ms sampled)\n"
        for call_key, call_entry in self.call_list.items():
            string += f"[{call_key}]: {call_entry.perf_time}ms ({call_entry.percent}%)\n"
        return string