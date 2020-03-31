import sys
sys.path.insert(0,'..')

from tqdm import trange
from itertools import product
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import os.path
import numpy as np
import pickle
import gillespy2
from gillespy2.solvers.python import *
#BasicSSASolver
from gillespy2.solvers.numpy import *
#BasicODESolver, BasicRootSolver, BasicTauLeapingSolver, NumPySSASolver, TauLeapingSolver
from gillespy2.solvers.cpp import *
#SSACSolver
from gillespy2.solvers.auto import *
#SSASolver
from gillespy2.solvers.stochkit import *
#StochKitODESolver, StochKitSolver


def timed_trials(models, solvers, trajectories, number_trials=30, override_number_trials={}, precompile_solvers=True,
                 output_file=None):
    """
    Runs a series of timed trials with a given set of GillesPy models and solvers.
    :param models: the list of models for the timed trials.
    :param solvers: the list of solvers for the timed trials.
    :param trajectories: the list of trajectories to time each model-solver combination on.
    :param number_trials: the number of times to run each model-solver-trajectory combination.
    :param override_number_trials: if specified, overrides the number of trials for specific solvers.
    :param precompile_solvers: whether or not to precompile solvers which may benefit from it (SSACSolver).
    :param output_file: if specified, file location where timing results will be stored using pickle.
    :return: A multi-dimensional dictionary which takes keys [model.name][solver.name] and returns NumPy arrays of the
    shape (len(trajectories), number_trials).
    """
    timing_data = {}
    for model in models:
        timing_data[model.name] = {}
    for model, solver in product(models, solvers):
        trials = number_trials
        if solver.name in override_number_trials:
            trials = override_number_trials[solver.name]
        if precompile_solvers and issubclass(solver, SSACSolver):
            solver = solver(model)
        timing_data[model.name][solver.name] = np.zeros((len(trajectories), 1+trials))
        for trajectory_i, trajectory in enumerate(trajectories):
            times = timing_data[model.name][solver.name][trajectory_i]
            times[0] = trajectory
            for i in trange(trials, desc = 'Model: {}, Solver: {}, Trajectories: {}'.format(model.name, solver.name, trajectory)):
                start = timer()
                model.run(solver=solver, number_of_trajectories=trajectory)
                stop = timer()
                times[1+i] = stop - start
            if output_file is not None:
                with open(output_file, 'wb') as f:
                    pickle.dump(timing_data, f)
    return timing_data


def plot_solver_run_times(timing_data, ylabel='Average seconds', reduce=np.mean, line_styles={}, transformation=None,
                          model_names=None, solver_names=None, baseline_solver_name=None, output_directory=None):
    """
    Plots matplotlib graphs comparing each solver's execution time for each model.
    :param timing_data: the timing data returned from a call to timed_trials().
    :param ylabel: label for the y axis.
    :param reduce: how to reduce the times for each model-solver-trajectory set, defaults to mean.
    :param line_styles: if specified, a dictionary denoting how to style lines for specific solvers.
    :param transformation: if specified, a function which will be called on the set of times for each solver.
    :param model_names: if specified, a list for a subset of the model names in the timing data to plot.
    :param solver_names: if specified, a list for a subset of the solver names in the timing data to plot.
    :param baseline_solver_name: if specified, timing results will be divided by the times of this solver.
    :param output_directory: if specified, plots will be saved in this folder.
    :return:
    """
    if model_names is None:
        model_names = timing_data.keys()
    for model in model_names:
        plt.figure()
        plt.title('{} Timing Results'.format(model))
        plt.xlabel('Trajectories')
        plt.ylabel(ylabel)
        baseline = None
        if baseline_solver_name is not None:
            baseline = reduce(timing_data[model][baseline_solver_name][:, 1:], axis=1)
            if transformation is not None:
                baseline = transformation(times)
        if solver_names is None:
            solver_names = timing_data[model].keys()
        for solver in solver_names:
            trajectories = timing_data[model][solver][:,0]
            times = reduce(timing_data[model][solver][:, 1:], axis=1)
            if transformation is not None:
                times = transformation(times)
            if baseline is not None:
                times = np.divide(baseline, times)
            if solver in line_styles:
                plt.plot(trajectories, times, line_styles[solver], label=solver)
            else:
                plt.plot(trajectories, times, label=solver)
        plt.legend(loc='best')
        if output_directory is not None:
            filename = '{}TimingResults{}.png'.format(model, ylabel.replace(' ', '_'))
            plt.savefig(os.path.join(output_directory, filename))
