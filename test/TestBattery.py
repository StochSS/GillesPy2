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

try:
    from tqdm import tqdm
except ImportError:
    raise ImportError('tqdm is required. Please install it.')
import statistics
from itertools import product
from scipy import stats
from timeit import default_timer as timer
import gillespy2
from gillespy2.solvers.numpy import *
# BasicODESolver, BasicRootSolver, BasicTauLeapingSolver, NumPySSASolver, TauLeapingSolver
from gillespy2.solvers.cpp import *
# SSACSolver
from gillespy2.solvers.auto import *
# SSASolver
from gillespy2.solvers.stochkit import *
# StochKitODESolver, StochKitSolver
from example_models import *


def timing_battery(number_of_samples, acceptable_deviation):
    gillespy2_home = ".."
    stochkit_home = "../../Stochkit/StochKit"
    output_file = "TestBattery_output.txt"

    solver_list = []
    key, value = None, None
    for key, value in globals().items():
        if isinstance(value, type) and issubclass(value, gillespy2.GillesPySolver) and value not in solver_list:
            solver_list.append(value)

    model_list = [create_decay(), create_trichloroethylene(), create_michaelis_menten(), create_schlogl()] #Update

    timing_list = []

    for model, solver in tqdm(product(model_list, solver_list)):
        median_list = []
        exterior_stats = []
        print("Testing Model : {} Solver : {}.".format(model.name, solver.name))
        for i in range(number_of_samples):
            try:
                solver = solver(model=model)
                standard_results = model.run(stochkit_home=stochkit_home)
                start = timer()
                test_results = model.run(solver=solver)
                stop = timer()
                median_list.append(stop-start)
                interior_stats = []
                for species in standard_results[0].keys():
                    if solver == NumPySSASolver:
                        deviation = stats.ks_2samp(standard_results[0][species], test_results[0][species])[0]
                        if deviation > acceptable_deviation:
                            print("Unacceptable deviation found on Model {} Solver {} with a deviation of {} on "
                                  "iteration {}. Exiting.".format(model.name, solver.name, deviation, i))
                            exit()
                        else:
                            interior_stats.append(deviation)
                    else:
                        deviation = stats.ks_2samp(standard_results[0][species], test_results[species])[0]
                        if deviation > acceptable_deviation:
                            print("Unacceptable deviation found on Model {} Solver {} with a deviation of {} on "
                                  "iteration {}. Exiting.".format(model.name, solver.name, deviation, i))
                            exit()
                        else:
                            interior_stats.append(deviation)
                exterior_stats.append(statistics.median(interior_stats))
            except Exception as e:
                print(e)
        median = statistics.median(median_list)
        timing_list.append([model, solver, median])
    with open(output_file, 'w') as out_file:
        out_file.write(*timing_list)


if __name__ == '__main__':
    timing_battery()

