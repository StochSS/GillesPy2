import sys
import json
import runpy
import shutil
import tempfile
import importlib

from pathlib import Path
from urllib import request
from zipfile import ZipFile
from collections import defaultdict
from argparse import ArgumentParser

from packaging import version
from packaging.version import Version

# Check to determine if a __init__.py file exists in the current test directory.
test_init = Path(__file__).parent.resolve().joinpath("test/__init__.py")
if not test_init.is_file():
    test_init.touch()

import gillespy2
import gillespy2.solvers.cpp
import gillespy2.solvers.numpy

from gillespy2.solvers.cpp import (
    ODECSolver,
    SSACSolver,
    TauLeapingCSolver
)
from gillespy2.solvers.numpy import (
    ODESolver,
    NumPySSASolver,
    TauLeapingSolver,
    TauHybridSolver
)

from test.example_models import *

target_models = [Tyson2StateOscillator, Example, VilarOscillator, MichaelisMenten, Dimerization]
profile_target_groups = {
    "c++": {
        ODECSolver:         target_models,
        SSACSolver:         target_models,
        TauLeapingCSolver:  target_models
    },
    "python": {
        ODESolver:          target_models,
        NumPySSASolver:     target_models,
        TauLeapingSolver:  target_models,
        TauHybridSolver:    target_models
    }
}

minimum_version = "1.6.0"
repository = "https://github.com/StochSS/GillesPy2"

def main():
    parser = ArgumentParser(description=f"The GillesPy2 per-version performance profiler. Versions > '{minimum_version}' are supported.")

    # --target-source and --target-version cannot both be set.
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument(
        "-t",
        "--target_source",
        action="store",
        dest="target_source",
        help=f"the directory path to a local copy of the target version repository. " 
                "This option takes precedence over '--target-version'. ",
        metavar="<target source>"
    )
    target_group.add_argument(
        "-v",
        "--target_version", 
        action="append",
        dest="target_versions",
        help=f"the GillesPy2 version to profile against. Version must be > '{minimum_version}'. "
                "Multiple versions can be specified with additional '-t VERSION' arguments.",
        metavar="<target version>"
    )

    # --group_name and --solver cannot both be set as they both specify solvers to profile.
    solver_group = parser.add_mutually_exclusive_group()
    solver_group.add_argument(
        "-g",
        "--group_name",
        action="store",
        dest="group_name",
        help="the name of the solver group to run. The possible options are 'c++', 'python', or 'all'. Defaults to 'all'.",
        metavar="<group name>",
        choices=[ "c++", "python", "all" ]
    )
    solver_group.add_argument(
        "-s",
        "--solver",
        action="append",
        dest="solver_names",
        help="solver name to profile. Multiple names can be specified with additional '-s SOLVER_NAME` arguments.",
        metavar="<solver name>",
        choices=[
            "ODECSolver",
            "SSACSolver",
            "TauLeapingCSolver",
            "ODESolver",
            "NumPySSASolver",
            "TauLeapingSolver",
            "TauHybridSolver"
        ],
    )

    parser.add_argument(
        "-n",
        "--number_of_runs",
        action="store",
        dest="run_count",
        help="the number of times each solver/model/version combo should be profiled.",
        metavar="<number of runs>",
        type=int,
        default=1
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output_file",
        help="the file profile results will be written to. If no value is passed then results will be written to stdout in JSON format.",
        metavar="<output file>",
        default=""
    )

    args = parser.parse_args()

    # If a value is not passed then these can be None.
    run_count = args.run_count
    output_file = args.output_file
    group_name = args.group_name
    solver_names = args.solver_names

    target_source = args.target_source
    target_versions = args.target_versions

    # Validate arguments.
    if run_count < 1:
        exit(
            f"!! Profile count of '{run_count}' is invalid because it is too small. "
            "Values must be greater than 0."
        )

    # Get the current GillesPy2 version.
    current_version = version.parse(gillespy2.__version__)
    targets = []

    # If the target_source directory was provided (and exists), grab the version number from it.
    if target_source is not None and Path(target_source).is_dir():
        # Grab the version of the specified source.
        target_source = Path(target_source).resolve()
        target_version = grab_local_version(target_source)

        # Validate the target_version against the current_version.
        validate_versions(current_version, target_version)

        targets.append((target_version, target_source))

    # If the target_source directory is not valid (or unset), clone the correct release by version.
    elif target_versions is not None and len(target_versions) > 0:
        for target_version in target_versions:
            # Validate the target_version against the current_version.
            target_version = version.parse(target_version)
            validate_versions(current_version, target_version)

            # Download the source of the specified version.
            target_source = clone_remote_version(target_version)

            # Get the local version of the downloaded source files.
            remote_version = grab_local_version(target_source)

            # Sanity check to ensure that the target and downloaded source files are of the same version.
            if target_version != grab_local_version(target_source):
                exit(
                    f"!! Failed to validate downloaded '{remote_version}' release files at '{target_source}' "
                    f"because the parsed local version is not the same as the target '{target_version}'."
                )

            # Save the target_version and target_source for later.
            targets.append((target_version, target_source))

    # If neither have valid inputs, quit.
    else:
        exit(
            "!! Failed to identify the target version as neither '--target-version' or '--target-source' are set. "
            "ONE of these must be set to a valid value to continue."
        )

    # Setup the profile targets.
    current = (current_version, Path(__file__).parent.resolve())
    profile_targets = { }

    # Determine the solver types to profile from the group_name value.
    if group_name is not None:
        print(f":: Using group '{group_name}'.")

        if group_name != "all":
            profile_targets = profile_target_groups[group_name.lower()]

    # If solver names were specified then determine the type and set to the target_models.
    elif solver_names is not None:
        for solver_name in solver_names:
            solver_type = getattr(gillespy2.solvers.cpp, solver_name) if hasattr(gillespy2.solvers.cpp, solver_name) else getattr(gillespy2.solvers.numpy, solver_name)

            profile_targets[solver_type] = target_models

    # If neither a group name or solver is specified set to default.
    else:
        profile_targets = { **profile_target_groups["c++"], **profile_target_groups["python"] }


    # Lets do this.
    print(f":: Ready to profile {current_version} against {', '.join(target_versions)}.")

    profile_results = run_profiler(current, targets, profile_targets, number_of_runs=run_count)
    results_json = json.dumps(profile_results, indent=4, default=vars)
    
    # If the output_file is not an empty string, write profiled results to it.
    if output_file != "":
        with Path(output_file).resolve().open("w") as outfile:
            outfile.write(results_json)

        return

    # Write the JSON results to stdout.
    print(results_json)


def validate_versions(current_version: Version, target_version: Version):
    # Check to ensure that the target version is not < the minimum version.
    if target_version < version.parse(minimum_version):
        exit(
            f"!! Failed to profile against target version '{target_version}' "
            f"because it is older than the minimum supported version '{minimum_version}'."
        )

def clone_remote_version(version: Version) -> Path:
    # Setup the download URL and working directory.
    download_url = f"{repository}/archive/refs/tags/v{version}.zip"

    temp_dir = tempfile.mkdtemp()
    destination = Path(temp_dir).joinpath("release.zip")
    source_dir = destination.parent.joinpath(f"GillesPy2-{version}")

    print(f":: Downloading version {version} from remote into '{destination}'...")

    # Download the release and stream it into the temp directory.
    dest_file = destination.open('wb')
    with request.urlopen(download_url) as incoming:
        shutil.copyfileobj(incoming, dest_file)

    # Verify that the downloaded file is where we expect it to be.
    if not destination.is_file():
        exit(
            f"!! Failed to download GillesPy2 release {version} to destination '{destination}'."
        )
    
    print(f" · :: Extracting '{destination}' into '{destination.parent}'...")

    # Extract the file into a subfolder of the temporary directory.
    with ZipFile(destination) as release:
        release.extractall(destination.parent)
    
    print(f" · :: GillesPy2 version {version} was successfully downloaded and extracted into '{source_dir}'.", end="\n\n")

    # This is a nasty patch to ensure that `test` directories are importable as modules.
    init_file = source_dir.joinpath("test", "__init__.py")
    if not init_file.is_file():
        init_file.touch()

    # The unzipped source files are now in the GillesPy2-version subdirectory.
    return source_dir

def grab_local_version(local_repo: Path) -> Version:
    version_file = local_repo.joinpath("gillespy2", "__version__.py")

    # Check that the version_file exists.
    if not version_file.is_file():
        exit(
            f"!! Failed to load version file at '{version_file}' "
            "because it does not exist. Please ensure that you are "
            "passing the correct path of the target source directory."
        )

    # Run the __version__.py file and save the results as a dict[str, Any].
    version_results = runpy.run_path(version_file)

    # Check that the __version__ property exists within the results of the prior call.
    if not "__version__" in version_results:
        exit(
            f"!! Failed to parse semantic version from version file at '{version_file}'. "
            "Please ensure that this file contains a valid '__version__' attribute.",
        )

    return version.parse(version_results["__version__"])

def run_profiler(
    current: "tuple[Version, Path]",
    targets: "list[tuple[Version, Path]]",
    profile_runs: "dict[type, list[type]]",
    number_of_runs: int = 1,
    ) -> "dict[str, dict[str, dict[str, list[object]]]]":

    # Store the results in a dictionary (for now).
    tree = lambda: defaultdict(tree)
    profile_results: "dict[str, dict[str, dict[str, list[object]]]]" = tree()

    # Determine profile targets. This makes iteration easier.
    profile_targets = [current] + targets

    # We need to wipe current imports so we can force in arbitrary GillesPy2 versions.
    for key in list(sys.modules.keys()):
        if key.startswith("gillespy2"):
            del sys.modules[key]

    for version, source in profile_targets:
        print(f":: Started profiler for version {version}.")

        # Enable modules to be loaded from the source directory.
        sys.path.insert(0, str(source))
        profiler = importlib.import_module("test.profiling.python_profiler", "test")

        # Disable GillesPy2 log output.
        core = importlib.import_module("gillespy2.core")
        core.log.disabled = True

        for solver, models in profile_runs.items():
            print("{:·<71} ::".format(f" · :: Profiling '{solver.__name__}' "))

            solver = importlib.import_module(solver.__module__).__getattribute__(solver.__name__)

            models = sorted(models * number_of_runs, key=lambda model: model.__name__)
            for model in models:
                message = f" ··· :: Running model '{model.__name__}' "
                print(message, end="", flush=True)

                if not isinstance(profile_results[solver.__name__][model.__name__][str(version)], list):
                    profile_results[solver.__name__][model.__name__][str(version)] = list()

                # Import required modules.
                model = importlib.import_module(model.__module__).__getattribute__(model.__name__)

                # Run the profiler and save the results.
                results = profiler.run_profiler(model=model(), solver=solver)
                profile_results[solver.__name__][model.__name__][str(version)].append(results)

                print("\r{:·<60} {:.2f} ms ::".format(message, results.execution_time))

            print("{:·<71} ::".format(f" · :: Done "), end="\n\n")

        # Remove the sys.path addition.
        sys.path.pop(0)

    return profile_results

if __name__ == "__main__":
    main()
