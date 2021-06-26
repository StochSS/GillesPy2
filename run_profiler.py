import sys
import runpy
import shutil
import tempfile
import importlib

from pathlib import Path
from urllib import request
from zipfile import ZipFile
from argparse import ArgumentParser

from packaging import version
from packaging.version import Version

# sys.path.append(Path(__file__).parent.joinpath("gillespy2").resolve())

import gillespy2 as gillespy2_current

minimum_version = "1.6.0"
repository = "https://github.com/StochSS/GillesPy2"

def main():
	parser = ArgumentParser(description=f"The GillesPy2 per-version performance profiler. Versions > '{minimum_version}' are supported.")
	parser.add_argument(
		"-s",
		"--target_source",
		action="store",
		dest="target_source",
		help=f"the directory path to a local copy of the target version repository. " 
		      "This option takes precedence over '--target-version'. "
	)
	parser.add_argument(
		"-t",
		"--target_version", 
		action="store",
		dest="target_version",
		help=f"the GillesPy2 version to profile against. Version must be > '{minimum_version}'."
	)

	args = parser.parse_args()

	# If a value is not passed then these can be None.
	target_source = args.target_source
	target_version = args.target_version

	# If the target_source directory was provided (and exists), grab the version number from it.
	if target_source is not None and Path(target_source).is_dir():
		target_source = Path(target_source).resolve()
		target_version = grab_local_version(target_source)

	# If the target_source directory is not valid (or unset), clone the correct release by version.
	elif target_version is not None:
		target_version = version.parse(target_version)
		target_source = clone_remote_version(target_version)

		# Get the local version of the downloaded source files.
		remote_version = grab_local_version(target_source)

		# Sanity check to ensure that the target and downloaded source files are of the same version.
		if target_version != grab_local_version(target_source):
			exit(
				f"Failed to validate downloaded '{remote_version}' release files at '{target_source}' "
				f"because the parsed local version is not the same as the target '{target_version}'."
			)

	# If neither have valid inputs, quit.
	else:
		exit(
			"Failed to identify the target version as neither '--target-version' or '--target-source' are set. "
			"ONE of these must be set to a valid value to continue."
		)

	# Get the current GillesPy2 version and validate.
	current_version = version.parse(gillespy2_current.__version__)
	validate_versions(current_version, target_version)

	# Lets do this.
	print(f"[INFO] Ready to profile target version '{target_version}' against current version '{current_version}'.")

	run_profiler(target_source)

def validate_versions(current_version: Version, target_version: Version):
	# Check to ensure that the target version is not < the minimum version.
	if target_version < version.parse(minimum_version):
		exit(
			f"Failed to profile against target version '{target_version}' "
			f"because it is older than the minimum supported version '{minimum_version}'."
		)

	# Get the current version of GillesPy2.
	current_version = version.parse(gillespy2_current.__version__)

	if current_version == target_version:
		exit(
			f"Failed to profile against target version '{target_version}' "
			"because it is identical to the current GillesPy2 version."
		)

def clone_remote_version(version: Version) -> Path:
	# Setup the download URL and working directory.
	download_url = f"{repository}/archive/refs/tags/v{version}.zip"

	temp_dir = tempfile.mkdtemp()
	destination = Path(temp_dir).joinpath("release.zip")
	source_dir = destination.parent.joinpath(f"GillesPy2-{version}")

	print(f"[INFO] Downloading version '{version}' from remote into '{destination}'...")

	# Download the release and stream it into the temp directory.
	dest_file = destination.open('wb')
	with request.urlopen(download_url) as incoming:
		shutil.copyfileobj(incoming, dest_file)

	# Verify that the downloaded file is where we expect it to be.
	if not destination.is_file():
		exit(
			f"Failed to download GillesPy2 release '{version}' to destination '{destination}'."
		)
	
	print(f"[INFO] Extracting '{destination}' into '{destination.parent}'...")

	# Extract the file into a subfolder of the temporary directory.
	with ZipFile(destination) as release:
		release.extractall(destination.parent)
	
	print(f"[INFO] GillesPy2 version '{version}' was successfully downloaded and extracted into '{source_dir}'.")

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
			f"Failed to load version file at '{version_file}' "
			"because it does not exist. Please ensure that you are "
			"passing the correct path of the target source directory."
		)

	# Run the __version__.py file and save the results as a dict[str, Any].
	version_results = runpy.run_path(version_file)

	# Check that the __version__ property exists within the results of the prior call.
	if not "__version__" in version_results:
		exit(
			f"Failed to parse semantic version from version file at '{version_file}'. "
			"Please ensure that this file contains a valid '__version__' attribute.",
		)

	return version.parse(version_results["__version__"])

def run_profiler(target_source: Path):
	# We need to generate a list of solvers and models to run our profilers with.
	profile_runs = {
		"NumPySSASolver": "Tyson2StateOscillator",
		"ODESolver": "Oregonator",
		"TauLeapingSolver": "VilarOscillator"
	}

	# We need to wipe current imports so we can force in arbitrary GillesPy2 versions.
	for key in list(sys.modules.keys()):
		if key.startswith("gillespy2"):
			del sys.modules[key]

	# Run the target version first.
	sys.path.insert(0, str(target_source))

	python_profiler = importlib.import_module("test.profiling.python_profiler", "test")
	example_models = importlib.import_module("test.example_models", "test")
	solvers = importlib.import_module("gillespy2.solvers.numpy", "test")

	print(
		"[DEBUG] Validating import paths:\n"
		f"\tpython_profiler:	'{python_profiler.__file__}'\n"
		f"\texample_models:		'{example_models.__file__}'\n"
		f"\tsolvers:		'{solvers.__file__}'\n"
	)

	for solver_name, model_name in profile_runs.items():
		model = getattr(example_models, model_name)()
		solver = getattr(solvers, solver_name)

		print(f"[INFO] Starting profile run with solver '{solver.name}' and model '{model.name}'...")
		results = python_profiler.run_profiler(model=model, solver=solver)

		print(results)

if __name__ == "__main__":
	main()
