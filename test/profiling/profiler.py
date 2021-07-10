import sys
import time
import copy
import runpy
import shutil
import tempfile

from pathlib import Path
from urllib import request
from zipfile import ZipFile
from collections import defaultdict
from packaging.version import Version

import plotly.graph_objects as plot_go

import gillespy2

from gillespy2.core import GillesPySolver, Model
from gillespy2.solvers import (
    ODECSolver,
    SSACSolver,
    TauLeapingCSolver
)

from test.example_models import *
from test.profiling.performance_data import PerformanceData

class Profiler:
    __target_models = [Tyson2StateOscillator, Example, VilarOscillator, MichaelisMenten, Dimerization]
    __default_runs = {
        ODECSolver:         __target_models,
        SSACSolver:         __target_models,
        TauLeapingCSolver:  __target_models
    }

    def __init__(self, targets: "list[ProfilerTarget]", runs: "dict[GillesPySolver, list[Model]]" = None):
        self.targets = copy.deepcopy(targets)
        self.runs = copy.deepcopy(runs)

        if self.runs is None:
            self.runs = self.__default_runs

        # Add the current version to the target list.
        self.current_version = gillespy2.__version__
        self.current_source = Path(gillespy2.__file__).parent.parent.resolve()

        # Sort the targets.
        self.targets = sorted(self.targets, key=lambda target: str(target.version))

    def run(self, number_of_runs: int = 1):
        tree = lambda: defaultdict(tree)
        results: "dict[GillesPySolver, dict[Model, dict[Version, list[PerformanceData]]]]" = tree()

        saved_path = copy.deepcopy(sys.path)

        # Wipe current gillespy2 imports.
        for key in list(sys.modules.keys()):
            if key.startswith("gillespy2") or key.startswith("test"):
                del sys.modules[key]

        for target in self.targets:
            version = target.version
            source = target.source

            print(f":: Started profiler for version {version}.")

            # Enable modules to be loaded from the source directory.
            sys.path.insert(0, str(source))

            import gillespy2.core
            gillespy2.core.log.disabled = True

            for solver, models in self.runs.items():
                print(f" · :: Profiling '{solver.__name__}'")

                # The type of solver to use needs to be detetmined based on the parent module of the solver.
                if "gillespy2.solvers.cpp" in solver.__module__:
                    from test.profiling import c_profiler as profiler
                
                else:
                    from test.profiling import python_profiler as profiler

                # solver = importlib.import_module(solver.__module__).__getattribute__(solver.__name__)

                models = sorted(models * number_of_runs, key=lambda model: model.__name__)
                for model in models:
                    message = f" ··· :: Running model '{model.__name__}' "
                    print(message, end="", flush=True)

                    if not isinstance(results[solver.__name__][model.__name__][str(version)], list):
                        results[solver.__name__][model.__name__][str(version)] = list()

                    # Import required modules.
                    # model = importlib.import_module(model.__module__).__getattribute__(model.__name__)

                    # Run the profiler and save the results.
                    start_time = time.perf_counter()
                    profile_results = profiler.run_profiler(model=model(), solver=solver)

                    stopwatch_time = (time.perf_counter() - start_time) * 1000
                    setattr(profile_results, "stopwatch_time", stopwatch_time)

                    # This weird assignment is needed since Python doesn't seem to recognize that `call_time` exists.
                    profile_results.call_list = profile_results.call_list

                    results[solver.__name__][model.__name__][str(version)].append(profile_results)

                    print("\r{:·<60} :: {:.2f} ms".format(message, profile_results.execution_time))

                print(f" · :: Done ", end="\n\n")

            # Remove the sys.path addition.
            sys.path.pop(0)

        # Restore sys.modules and sys.path.
        for key in list(sys.modules.keys()):
            if key.startswith("gillespy2") or key.startswith("test"):
                del sys.modules[key]

        sys.path = saved_path

        # Instantiate and return a ProfileResults object.
        return ProfilerResults(results, self.targets, self.runs, number_of_runs)

class ProfilerTarget:
    @classmethod
    def from_current(cls) -> "ProfilerTarget":
        local_path = Path(gillespy2.__file__).parent.parent.resolve()
        print(local_path)

        return cls.from_source(gillespy2.__version__, local_path)

    @classmethod
    def from_version(cls, version_str: str) -> "ProfilerTarget":
        version = Version(version_str)
        profiler_target = cls(version)

        local_source = profiler_target.download_source()

        return cls.from_source(version_str, local_source)

    @classmethod
    def from_source(cls, version_str: str, local_source: Path) -> "ProfilerTarget":
        version = Version(version_str)
        profiler_target = cls(version)

        # Validate the source.
        profiler_target.validate_source(local_source)
        profiler_target.source = local_source

        return profiler_target

    def __init__(self, version: Version):
        self.minimum_version = Version("1.6.0")
        self.current_version = Version(gillespy2.__version__)

        self.version = version
        self.source = None

        if self.version < self.minimum_version:
            raise Exception(
                f"Failed to insantiate target with version '{self.version}' "
                f"because it is older than the minimum supported version '{self.minimum_version}'."
            )

        if self.version > self.current_version:
            raise Exception(
                f"Failed to instantiate target with version '{self.version}' "
                f"because it is newer than the current version '{self.current_version}'."
            )
    
    def validate_source(self, source_dir: Path):
        version_file = source_dir.joinpath("gillespy2", "__version__.py")

        # The version file needs to exist.
        if not version_file.is_file():
            raise Exception( 
                f"Failed to load version file at '{version_file}' "
                "because it does not exist. Please ensure that you are "
                "passing the correct path to the target directory."
            )

        # Run the __version__.py file and save the results.
        version_results = runpy.run_path(str(version_file))

        # Check that the __version__property exists within the results.
        if not "__version__" in version_results:
            raise Exception(
                f"Failed to parse semantic version from version file at '{version_file}'. "
                "Please ensure that this file contains a valid '__version__' attribute."
            )

        if Version(version_results["__version__"]) != self.version:
            raise Exception(
                "Failed while validating source version -- the target and source versions do not match."
            )

        test_init = Path(source_dir).joinpath("test/__init__.py")
        if not test_init.is_file():
            test_init.touch()

        # If the version is 1.6.0 to 1.6.2 then an updated Makefile needs to be hot-patched in. These versions
        # use a Makefile inside of the 'profiling' directory which does not include build instructions for the 
        # ArgParser.
        if self.version >= Version("1.6.0") and self.version <= Version("1.6.2"):
            # The path of the profiling Makefile.
            old_makefile = source_dir.joinpath("test/profiling/Makefile")

            if not old_makefile.is_file():
                print(f"Expected to find old Makefile, found nothing. This is undefined behavior: {old_makefile}.")

            # Remove the old makefile and replace with patched version.
            old_makefile.unlink()
            new_makefile = Path(__file__).parent.joinpath("assets/1_6_0-1_6_2-makefile")

            shutil.copyfile(new_makefile, old_makefile)

            # We've patched the Makefile, but the c_profiler.py script still uses '-' prefixed arguments.
            # This is a quick and nasty patch.
            profiler_path = Path(source_dir).joinpath("test/profiling/c_profiler.py")

            if not profiler_path.is_file():
                print(f"Expected to find c_profiler.py, found nothing. This is undefined behavior: {profiler_path}.")

            # Read the contents of the profiler_path file.
            with profiler_path.open("r") as infile:
                profiler_contents = infile.read()

            # Replace the following instances.
            needs_patching = ["-trajectories", "-timesteps", "-end", "-increment"]

            for target in needs_patching:
                profiler_contents = profiler_contents.replace(target, target.replace("-", "--"))

            # Write the patched contents back to the file.
            profiler_path.unlink()

            with profiler_path.open("w") as outfile:
                outfile.write(profiler_contents)

            print(f"Successfully patched {self.version}.")

    def download_source(self) -> Path:
        download_url = f"https://github.com/StochSS/GillesPy2/archive/refs/tags/v{self.version}.zip"

        temp_dir = tempfile.mkdtemp()
        destination = Path(temp_dir).joinpath("release.zip")
        source_dir = destination.parent.joinpath(f"GillesPy2-{self.version}")

        # Download the release and stream it into the temp directory.
        dest_file = destination.open('wb')
        with request.urlopen(download_url) as incoming:
            shutil.copyfileobj(incoming, dest_file)

        # Verify that the downloaded file is where we expect it to be.
        if not destination.is_file():
            raise Exception(
                f"Failed to download GillesPy2 release {self.version} to destination '{destination}'."
            )

        # Extract the file into a subfolder of the temporary directory.
        with ZipFile(destination) as release:
            release.extractall(destination.parent)

        # This is a nasty patch to ensure that `test` directories are importable as modules.
        init_file = source_dir.joinpath("test", "__init__.py")
        if not init_file.is_file():
            init_file.touch()

        # The unzipped source files are now in the GillesPy2-version subdirectory.
        return source_dir

class ProfilerResults:
    def __init__(self, results, targets, runs, run_count):
        self.line_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]

        self.results: "dict[GillesPySolver, dict[Model, dict[Version, list[PerformanceData]]]]" = results
        self.targets = targets
        self.runs = runs
        self.run_count = run_count

    def plot_solver(self, solver: "GillesPySolver", target_stat: str = "execution_time"):
        # Generate lists of model and version names.
        models = [(model.__name__) for model in self.runs[solver]]
        versions = [(str(target.version)) for target in self.targets]

        bars = []
        for index, version in enumerate(versions):
            times = []

            for model in models:
                all_times = [(getattr(profile_data, target_stat, 0.0)) for profile_data in self.results[solver.__name__][model][version]]
                average = round(sum(all_times) / len(all_times), 2)

                times.append(average)

            bars.append(plot_go.Bar(
                name=f"Version {version}", 
                x=models, 
                y=times, 
                text=times, 
                textposition="auto", 
                marker_color=self.line_colors[index]))

        figure = plot_go.Figure(data=bars)
        figure.update_layout(
            barmode="group", 
            title=f"{solver.__name__}",
            xaxis_tickfont_size=14,
            yaxis=dict(
                title=f"{target_stat.replace('_', ' ').capitalize()} (ms)",
                titlefont_size=16,
                tickfont_size=14
            ),
            legend=dict(
                x=0,
                y=1.0,
                bgcolor="rgba(255, 255, 255, 0)",
                bordercolor="rgba(255, 255, 255, 0)"
            ))

        figure.show()
