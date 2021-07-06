import sys
import copy
import runpy
import shutil
import tempfile
import importlib

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

        # If a profile target for the current version already exists, remove it.
        self.targets = [(target) for target in self.targets if target.version != self.current_version]
        self.targets.insert(0, ProfilerTarget.from_source(self.current_version, self.current_source))

        # Sort the targets.
        self.targets = sorted(self.targets, key=lambda target: str(target.version))

    def run(self, number_of_runs: int = 1):
        tree = lambda: defaultdict(tree)
        results: "dict[GillesPySolver, dict[Model, dict[Version, list[PerformanceData]]]]" = tree()

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

                solver = importlib.import_module(solver.__module__).__getattribute__(solver.__name__)

                models = sorted(models * number_of_runs, key=lambda model: model.__name__)
                for model in models:
                    message = f" ··· :: Running model '{model.__name__}' "
                    print(message, end="", flush=True)

                    if not isinstance(results[solver.__name__][model.__name__][str(version)], list):
                        results[solver.__name__][model.__name__][str(version)] = list()

                    # Import required modules.
                    model = importlib.import_module(model.__module__).__getattribute__(model.__name__)

                    # Run the profiler and save the results.
                    profile_results = profiler.run_profiler(model=model(), solver=solver)

                    # This weird assignment is needed since Python doesn't seem to recognize that `call_time` exists.
                    profile_results.call_list = profile_results.call_list
                    results[solver.__name__][model.__name__][str(version)].append(profile_results)

                    print("\r{:·<60} :: {:.2f} ms".format(message, profile_results.execution_time))

                print(f" · :: Done ", end="\n\n")

            # Remove the sys.path addition.
            sys.path.pop(0)

        # Instantiate and return a ProfileResults object.
        return ProfilerResults(results, self.targets, self.runs, number_of_runs)

class ProfilerTarget:
    @classmethod
    def from_source(cls, version: str, local_source: Path = None) -> "ProfilerTarget":
        version = Version(version)
        profiler_target = cls(version)
        
        # If the local_source is None, download the repository.
        if local_source is None:
            local_source = profiler_target.download_source()

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
        version_results = runpy.run_path(version_file)

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

    def plot_solver(self, solver: "GillesPySolver"):
        # Generate lists of model and version names.
        models = [(model.__name__) for model in self.runs[solver]]
        versions = [(str(target.version)) for target in self.targets]

        bars = []
        for index, version in enumerate(versions):
            times = []

            for model in models:
                all_times = [(profile_data.execution_time) for profile_data in self.results[solver.__name__][model][version]]
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
                title="Execution time (ms)",
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
