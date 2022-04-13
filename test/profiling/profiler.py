import argparse

import c_profiler
import performance_data
import gillespy2
from examples.StartingModels.VilarOscillator.VilarOscillator import VilarOscillator


def run_profiler(model: "gillespy2.Model", solver_name: "str",
                 iterations=1,
                 include_run_times=False,
                 include_call_times=False):
    # Validate/override call args
    if not isinstance(iterations, int) or iterations <= 0:
        iterations = 1
    
    ensemble = performance_data.PerformanceEnsemble()

    print(f"=== === ===  [Solver = {solver_name}, Model = {model.name}]  === === ===")
    for run_number, perf_results in enumerate(
            c_profiler.run_profiler(model=model, solver=solver_name, iterations=iterations)):
        print(f"[Run {run_number+1} out of {iterations}]")
        ensemble.runs.append(perf_results)
        if not include_run_times:
            continue

        print(f" Execution Time:\t{perf_results.execution_time:0.1f}ms")
        print(f" CPU Time:\t\t{perf_results.sample_time:0.1f}ms")
        if include_call_times:
            for entry, time in list(perf_results.call_list.items())[:10]:
                print(f" * {time.perf_time:0.1f}ms:\t{entry}")

    avg_time = ensemble.get_average_time()
    total_time = ensemble.get_total_time()
    print(f"[Ensemble Results]")
    print(f" Average Execution Time: {avg_time.execution_time:0.2f}ms ({total_time.execution_time:0.2f}ms Total)")
    print(f" Average CPU Time:       {avg_time.sample_time:0.2f}ms ({total_time.sample_time:0.2f}ms Total)")


if __name__ == "__main__":
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument("-n", "--runs", type=int, default=1, help="number of times to run the profiler")
        parser.add_argument("-c", "--call-times", action="store_true", help="include call times for particular function calls")
        parser.add_argument("--ensemble-only", action="store_true", help="do not include data for each run, only ensemble data")
        return parser.parse_args()

    args = parse_args()
    run_profiler(
        model=VilarOscillator(),
        solver_name=gillespy2.SSACSolver.target,
        include_run_times=(not args.ensemble_only),
        iterations=args.runs,
        include_call_times=args.call_times,
    )
