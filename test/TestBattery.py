from tqdm import tqdm
import click
import statistics
from itertools import product
from scipy import stats
from timeit import default_timer as timer


@click.command()
@click.option('--gillespy2_home', default="/home/jackson/Research/GillesPy2/", help='Location of Gillespy2 Directory')
@click.option('--number_of_samples', prompt='Sample Size', help='The Number of times that you want to run it.')
@click.option('--stochkit_home', prompt='Location of Stochkit', help='The Number of times that ')
@click.option('--acceptable deviation', prompt='Numeric Value of Statistical Variation', help='')
@click.option('--output_file', prompt='Name of Output File', help='The name of the output file')
def timing_battery(gillespy2_home, number_of_samples, stochkit_home, acceptable_deviation, output_file):


    from .basic_ssa_solver import BasicSSASolver
    from .basic_root_solver import BasicRootSolver
    from .optimized_ssa_solver import OptimizedSSASolver

    from Repository.gillespy2.example_models import *

    model_list = [Example(), Trichloroethylene(), MichaelisMenten(), Schlogl()]

    solver_list = [BasicSSASolver, BasicRootSolver, OptimizedSSASolver, SSACSolver()]
    timing_list = []

    for model, solver in tqdm(product(model_list, solver_list)):
        median_list = []
        exterior_stats = []
        print("Testing Model : {} Solver : {}.".format(model.name, solver.name))
        for i in range(number_of_samples):
            try:
                standard_results = model.run(stochkit_home=stochkit_home)
                start = timer()
                test_results = model.run(solver=solver)
                stop = timer()
                median_list.append(stop-start)
                interior_stats = []
                for species in standard_results[0].keys():
                    if solver == OptimizedSSASolver:
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

