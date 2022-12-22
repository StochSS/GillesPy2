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
import os
import csv

from datetime import datetime
from collections import UserDict, UserList

import numpy as np

from gillespy2.core.jsonify import Jsonify
from gillespy2.core.gillespyError import ValidationError

def common_rgb_values():
    '''
    List of 50 hex color values used for plotting graphs
    '''
    return [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
        '#bcbd22', '#17becf', '#ff0000', '#00ff00', '#0000ff', '#ffff00', '#00ffff', '#ff00ff',
        '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#ff9999', '#ffcc99',
        '#ccff99', '#cc99ff', '#ffccff', '#62666a', '#8896bb', '#77a096', '#9d5a6c', '#9d5a6c',
        '#eabc75', '#ff9600', '#885300', '#9172ad', '#a1b9c4', '#18749b', '#dadecf', '#c5b8a8',
        '#000117', '#13a8fe', '#cf0060', '#04354b', '#0297a0', '#037665', '#eed284', '#442244',
        '#ffddee', '#702afb'
    ]

def _plot_iterate(self, show_labels=True, included_species_list=[]):
    import matplotlib.pyplot as plt # pylint: disable=import-outside-toplevel
    for i, species in enumerate(self.data):
        if species != 'time':

            if species not in included_species_list and included_species_list:
                continue

            line_color = common_rgb_values()[(i - 1) % len(common_rgb_values())]

            if show_labels:
                label = species
            else:
                label = ""

            plt.plot(self.data['time'], self.data[species], label=label, color=line_color)


def _plotplotly_iterate(trajectory, show_labels=True, trace_list=None, line_dict=None, included_species_list=[]):
    """
    Helper method for Results .plotplotly() method
    """

    if trace_list is None:
        trace_list = []

    import plotly.graph_objs as go # pylint: disable=import-outside-toplevel

    for i, species in enumerate(trajectory.data):
        if species != 'time':

            if species not in included_species_list and included_species_list:
                continue

            if line_dict is None:
                line_dict = {}

            # If number of species exceeds number of available colors, loop back through colors
            line_dict['color'] = common_rgb_values()[(i-1)%len(common_rgb_values())]

            if show_labels:
                trace_list.append(
                    go.Scatter(
                        x=trajectory.data['time'],
                        y=trajectory.data[species],
                        mode='lines',
                        name=species,
                        line=line_dict,
                        legendgroup=species
                    )
                )
            else:
                trace_list.append(
                    go.Scatter(
                        x=trajectory.data['time'],
                        y=trajectory.data[species],
                        mode='lines',
                        name=species,
                        line=line_dict,
                        legendgroup=species,
                        showlegend=False
                    )
                )

    return trace_list

class Trajectory(UserDict, Jsonify):
    """ Trajectory Dict created by a gillespy2 solver containing single trajectory, extends the UserDict object.

    :param data: A dictionary of trajectory values created by a solver
    :type data: UserDict

    :param model: The name of the model used to create the trajectory
    :type model: str

    :param solver_name: The name of the solver used to create the trajectory
    :type solver_name: str

    :param rc: The solvers status return code.
    :type rc: int

    :param status: The solver status ('Success','Timed out')
    """

    def __init__(self, data, model=None, solver_name="Undefined solver name", rc=0):

        self.data = data
        self.model = model
        self.solver_name = solver_name
        self.rc = rc

        status_list = {0: 'Success', 33: 'Timed Out'}
        self.status = status_list[rc]

    def __getitem__(self, key):
        if isinstance(key, int):
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            species = list(self.data.keys())[key]
            msg = "Trajectory is of type dictionary."
            msg += f"Use trajectory['[{species}]'] instead of trajectory[{key}]['{species}']"
            msg += f"Retrieving trajectory['[{species}]']"
            log.warning(msg)
            return self.data[species]
        if key in self.data:
            return self.data[key]
        if hasattr(self.__class__, "__missing__"):
            return self.__class__.__missing__(self, key)
        raise KeyError(key)


class Results(UserList, Jsonify):
    """
    List of Trajectory objects created by a gillespy2 solver, extends the UserList object.

    :param data: A list of trajectory objects
    :type data: UserList
    """

    def __init__(self, data):
        self.data = data

    def __getattribute__(self, key):
        if key in ('model', 'solver_name', 'rc', 'status'):
            if len(self.data) > 1:
                from gillespy2.core import log # pylint: disable=import-outside-toplevel
                msg = f"Results is of type list. Use results[i]['{key}'] instead of results['{key}']"
                log.warning(msg)
            return getattr(Results.__getattribute__(self, key='data')[0], key)
        return UserList.__getattribute__(self, key)

    def __getitem__(self, key):
        if key == 'data':
            return UserList.__getitem__(self, key)
        if isinstance(key, str):
            if len(self.data) > 1:
                from gillespy2.core import log # pylint: disable=import-outside-toplevel
                msg = f"Results is of type list. Use results[i]['{key}'] instead of results['{key}']"
                log.warning(msg)
            return self.data[0][key]
        return UserList.__getitem__(self,key)

    def __add__(self, other):
        combined_data = Results(data=(self.data + other.data))
        consistent_solver = combined_data._validate_solver()
        consistent_model = combined_data._validate_model()

        if not consistent_solver:
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            log.warning("Results objects contain Trajectory objects from multiple solvers.")

        if not consistent_model:
            raise ValidationError('Results objects contain Trajectory objects from multiple models.')

        return combined_data

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other)

    def _validate_model(self, reference=None):
        is_valid = True
        if reference is not None:
            reference_model = reference
        else:
            reference_model = self.data[0].model.get_json_hash()
        for trajectory in self.data:
            if trajectory.model.get_json_hash() != reference_model:
                is_valid = False
        return is_valid

    def _validate_solver(self, reference=None):
        is_valid = True
        if reference is not None:
            reference_solver = reference
        else:
            reference_solver = self.data[0].solver_name
        for trajectory in self.data:
            if trajectory.solver_name != reference_solver:
                is_valid = False
        return is_valid

    def _validate_title(self, show_title):
        if not show_title:
            title = ''
            return title
        if self._validate_model():
            title_model = self.data[0].model.name
        else:
            title_model = 'Multiple Models'
        if self._validate_solver():
            title_solver = self.data[0].solver_name
        else:
            title_solver = 'Multiple Solvers'
        title = (title_model + " - " + title_solver)
        return title

    def to_array(self):
        '''
        Convert the results object into a numpy array.

        :returns: Array containing the result of the simulation.
        :rtype: numpy.ndarray
        '''
        results = []
        size1 = len(self.data[0]['time'])
        size2 = len(self.data[0])

        for trajectory in self.data:
            new_array = np.zeros((size1, size2))
            for j, key in enumerate(trajectory):
                new_array[:, j] = trajectory[key]
            results.append(new_array)
        return np.array(results)

    @classmethod
    def build_from_solver_results(cls, solver, live_output_options):
        """
        Build a gillespy2.Results object using the provided solver results.

        :param solver: The solver used to run the simulation.
        :type solver: gillespy2.GillesPySolver

        :param live_output_options: dictionary contains options for live_output. By default {"interval":1}.
            "interval" specifies seconds between displaying.
            "clear_output" specifies if display should be refreshed with each display
        :type live_output_options: dict
        """
        if solver.rc == 33:
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            log.warning('GillesPy2 simulation exceeded timeout.')
        if hasattr(solver.result[0], 'shape'):
            return solver.result
        if len(solver.result) > 0:
            results_list = []
            for result in solver.result:
                temp = Trajectory(data=result, model=solver.model, solver_name=solver.name, rc=solver.rc)
                results_list.append(temp)

            results = Results(results_list)
            if "type" in live_output_options.keys() and live_output_options['type'] == "graph":
                results.plot()
            return results
        raise ValueError("number_of_trajectories must be non-negative and non-zero")

    def to_csv(self, path=None, nametag=None, stamp=None):
        """
        Outputs the Results to one or more .csv files in a new directory.

        :param nametag: allows the user to optionally "tag" the directory and included files. Defaults to the model
            name.
        :type nametag: str

        :param path: the location for the new directory and included files. Defaults to model location.
        :type path: str

        :param stamp: Allows the user to optionally "tag" the directory (not included files). Default is timestamp.
        :type stamp: str
        """
        if stamp is None:
            now = datetime.now()
            stamp = datetime.timestamp(now)
        if nametag is None:
            identifier = self._validate_title(show_title=True)
        else:
            identifier = nametag
        if path is None:
            directory = os.path.join(".", str(identifier)+str(stamp))
        else:
            directory = os.path.join(path, str(identifier)+str(stamp))
        # multiple trajectories
        if isinstance(self.data, list):
            os.mkdir(directory)
            for i, trajectory in enumerate(self.data):  # write each CSV file
                filename = os.path.join(directory, str(identifier)+str(i)+".csv")
                field_names = []
                for species in trajectory:  # build the header
                    field_names.append(species)
                with open(filename, 'w', newline='', encoding="utf-8") as csv_file:
                    csv_writer = csv.writer(csv_file)
                    csv_writer.writerow(field_names)  # write the header
                    for j, _ in enumerate(trajectory['time']):  # write all lines of the CSV file
                        this_line=[]
                        for species in trajectory:  # build one line of the CSV file
                            this_line.append(trajectory[species][j])
                        csv_writer.writerow(this_line)  # write one line of the CSV file

    def plot(self, index=None, xaxis_label="Time", xscale='linear', yscale='linear', yaxis_label="Value",
             style="default", title=None, show_title=False, show_legend=True, multiple_graphs=False,
             included_species_list=[], save_png=False, figsize=(18, 10)):
        """
        Plots the Results using matplotlib.

        :param index: If not none, the index of the Trajectory to be plotted.
        :type index: int

        :param xaxis_label: The label for the x-axis
        :type xaxis_label: str

        :param yaxis_label: The label for the y-axis
        :type yaxis_label: str

        :param title: The title of the graph
        :type title: str

        :param multiple_graphs: IF each trajectory should have its own graph or if they should overlap.
        :type multiple_graphs: bool

        :param included_species_list: A list of strings describing which species to include. By default displays all
            species.
        :type included_species_list: list

        :param save_png: Should the graph be saved as a png file. If True, File name is title of graph. If a string is
            given, file is named after that string.
        :type save_png: bool or str

        :param figsize: The size of the graph. A tuple of the form (width,height). Is (18,10) by default.
        :type figsize: tuple of ints (x,y)
        """

        import matplotlib.pyplot as plt # pylint: disable=import-outside-toplevel
        from collections.abc import Iterable # pylint: disable=import-outside-toplevel
        trajectory_list = []
        if isinstance(index, Iterable):
            for i in index:
                trajectory_list.append(self.data[i])
        elif isinstance(index, int):
            trajectory_list.append(self.data[index])
        else:
            trajectory_list = self.data

        if title is None:
            title = self._validate_title(show_title)

        if len(trajectory_list) < 2:
            multiple_graphs = False

        if multiple_graphs:
            for i, trajectory in enumerate(trajectory_list):
                result = Results(data=[trajectory])
                if isinstance(save_png, str):
                    result.plot(xaxis_label=xaxis_label, yaxis_label=yaxis_label, title=title + " " + str(i + 1),
                                style=style, included_species_list=included_species_list, save_png=save_png + str(i + 1)
                                , figsize=figsize)
                else:
                    result.plot(xaxis_label=xaxis_label, yaxis_label=yaxis_label, title=title + " " + str(i + 1),
                                style=style, included_species_list=included_species_list, save_png=save_png,
                                figsize=figsize)

        else:
            try:
                plt.style.use(style)
            except Exception:
                from gillespy2.core import log # pylint: disable=import-outside-toplevel
                msg = f"Invalid matplotlib style. Try using one of the following {plt.style.available}"
                log.warning(msg)
                plt.style.use("default")

            plt.figure(figsize=figsize)
            plt.title(title, fontsize=18)
            plt.xlabel(xaxis_label)
            plt.ylabel(yaxis_label)
            plt.xscale(xscale)
            plt.yscale(yscale)

            for i, trajectory in enumerate(trajectory_list):

                if i > 0:
                    _plot_iterate(trajectory, included_species_list=included_species_list, show_labels=False)
                else:
                    _plot_iterate(trajectory, included_species_list=included_species_list)

            if show_legend:
                plt.legend(loc='best')
            plt.plot()
            if isinstance(save_png, str):
                plt.savefig(save_png)

            elif save_png:
                plt.savefig(title)

    def plotplotly(self, index=None, xaxis_label="Time", yaxis_label="Value", title=None,
                   show_title=False, show_legend=True, multiple_graphs=False, included_species_list=[],
                   return_plotly_figure=False,  **layout_args):
        """ Plots the Results using plotly. Can only be viewed in a Jupyter Notebook.

        :param index: If not none, the index of the Trajectory to be plotted.
        :type index: int

        :param xaxis_label: The label for the x-axis
        :type xaxis_label: str

        :param yaxis_label: The label for the y-axis
        :type yaxis_label: str

        :param title: The title of the graph
        :type title: str

        :param show_title: If True, title will be shown on graph.
        :type show_title: bool

        :param show_legend: Default True, if False, legend will not be shown on graph.
        :type show_legend: bool

        :param multiple_graphs: IF each trajectory should have its own graph or if they should overlap.
        :type multiple_graphs: bool

        :param included_species_list: A list of strings describing which species to include. By default displays all
            species.
        :type included_species_list: list

        :param return_plotly_figure: Whether or not to return a figure dictionary of data(graph object traces) and
            layout which may be edited by the user
        :type return_plotly_figure: bool

        :param **layout_args: Optional additional arguments to be passed to plotlys layout constructor.
        :type **layout_args: dict
        """

        from plotly.offline import init_notebook_mode, iplot # pylint: disable=import-outside-toplevel
        import plotly.graph_objs as go # pylint: disable=import-outside-toplevel

        # Backwards compatibility with xaxis_label argument (which duplicates plotly's xaxis_title argument)
        if layout_args.get('xaxis_title') is not None:
            xaxis_label = layout_args.get('xaxis_title')
            layout_args.pop('xaxis_title')
        if layout_args.get('yaxis_title') is not None:
            yaxis_label = layout_args.get('yaxis_title')
            layout_args.pop('yaxis_title')

        init_notebook_mode(connected=True)

        from collections.abc import Iterable # pylint: disable=import-outside-toplevel
        trajectory_list = []
        if isinstance(index, Iterable):
            for i in index:
                trajectory_list.append(self.data[i])
        elif isinstance(index, int):
            trajectory_list.append(self.data[index])
        else:
            trajectory_list = self.data

        number_of_trajectories = len(trajectory_list)

        if title is None:
            title = self._validate_title(show_title)

        fig = dict(data=[], layout=[])

        if len(trajectory_list) < 2:
            multiple_graphs = False

        if multiple_graphs:

            from plotly import subplots # pylint: disable=import-outside-toplevel

            fig = subplots.make_subplots(print_grid=False, rows=int(number_of_trajectories/2) +
                                                             int(number_of_trajectories % 2), cols=2)

            for i, trajectory in enumerate(trajectory_list):
                if i > 0:
                    trace_list = _plotplotly_iterate(trajectory, trace_list=[], included_species_list=
                    included_species_list, show_labels=False)
                else:
                    trace_list = _plotplotly_iterate(trajectory, trace_list=[], included_species_list=
                    included_species_list)

                for trace in trace_list:
                    if i % 2 == 0:
                        fig.append_trace(trace, int(i/2) + 1, 1)
                    else:
                        fig.append_trace(trace, int(i/2) + 1, 2)

                fig['layout'].update(
                    autosize=True, height=400*len(trajectory_list), showlegend=show_legend, title=title
                )
        else:
            trace_list = []
            for i, trajectory in enumerate(trajectory_list):
                if i > 0:
                    trace_list = _plotplotly_iterate(trajectory, trace_list=trace_list, included_species_list=
                    included_species_list, show_labels=False)
                else:
                    trace_list = _plotplotly_iterate(trajectory, trace_list=trace_list, included_species_list=
                    included_species_list)


            layout = go.Layout(
                showlegend=show_legend,
                title=title,
                xaxis_title=xaxis_label,
                yaxis_title=yaxis_label,
                **layout_args
            )

            fig['data'] = trace_list
            fig['layout'] = layout

        if not return_plotly_figure:
            iplot(fig)
            return None
        return fig

    def average_ensemble(self):
        """
        Generate a single Results object with a Trajectory that is made of the means of all trajectories' outputs

        :returns: The Results object
        """

        trajectory_list = self.data
        number_of_trajectories = len(trajectory_list)

        output_trajectory = Trajectory(data={}, model=trajectory_list[0].model, solver_name=
        trajectory_list[0].solver_name)

        for species in trajectory_list[0]:  # Initialize the output to be the same size as the inputs
            output_trajectory[species] = [0]*len(trajectory_list[0][species])

        output_trajectory['time'] = trajectory_list[0]['time']

        # Add every value of every Trajectory Dict into one output Trajectory
        for i in range(0, number_of_trajectories):
            trajectory_dict = trajectory_list[i]
            for species in trajectory_dict:
                if species == 'time':
                    continue
                for k in range(0, len(output_trajectory[species])):
                    output_trajectory[species][k] += trajectory_dict[species][k]

        for species in output_trajectory:   # Divide for mean of every value in output Trajectory
            if species == 'time':
                continue
            for i in range(0, len(output_trajectory[species])):
                output_trajectory[species][i] /= number_of_trajectories

        output_results = Results(data=[output_trajectory])  # package output_trajectory in a Results object

        return output_results

    def stddev_ensemble(self, ddof=0):
        """
        Generate a single Results object with a Trajectory that is made of the sample standard deviations of all
        trajectories' outputs.

        :param ddof: Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents
            the number of trajectories. Sample standard deviation uses ddof of 1. Defaults to population standard
            deviation where ddof is 0.
        :type ddof: int

        :returns: the Results object
        """

        from math import sqrt # pylint: disable=import-outside-toplevel

        trajectory_list = self.data
        number_of_trajectories = len(trajectory_list)

        if ddof == number_of_trajectories:
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            log.warning("ddof must be less than the number of trajectories. Using ddof of 0")
            ddof = 0

        average_list = self.average_ensemble().data[0]

        output_trajectory = Trajectory(data={}, model=trajectory_list[0].model, solver_name=
        trajectory_list[0].solver_name)

        for species in trajectory_list[0]:  # Initialize the output to be the same size as the inputs
            output_trajectory[species] = [0]*len(trajectory_list[0][species])

        output_trajectory['time'] = trajectory_list[0]['time']

        for i in range(0, number_of_trajectories):
            trajectory_dict = trajectory_list[i]
            for species in trajectory_dict:
                if species == 'time':
                    continue
                for k in range(0, len(output_trajectory['time'])):
                    output_trajectory[species][k] += (trajectory_dict[species][k] - average_list[species][k])\
                                          * (trajectory_dict[species][k] - average_list[species][k])

        for species in output_trajectory:   # Divide for mean of every value in output Trajectory
            if species == 'time':
                continue
            for i in range(0, len(output_trajectory[species])):
                output_trajectory[species][i] /= (number_of_trajectories - ddof)
                output_trajectory[species][i] = sqrt(output_trajectory[species][i])

        output_results = Results(data=[output_trajectory])  # package output_trajectory in a Results object
        return output_results


    def plotplotly_mean_stdev(self, xaxis_label="Time", yaxis_label="Value", title=None,
                                 show_title=False, show_legend=True, included_species_list=[],
                                 return_plotly_figure=False, ddof=0, **layout_args):
        """
        Plot a plotly graph depicting the mean and standard deviation of a results object

        :param xaxis_label: The label for the x-axis
        :type xaxis_label: str

        :param yaxis_label: The label for the y-axis
        :type yaxis_label: str

        :param title: The title of the graph
        :type title: str

        :param show_title: If True, title will be shown on graph.
        :type show_title: bool

        :param show_legend: Default True, if False, legend will not be shown on graph.
        :type show_legend: bool

        :param included_species_list: A list of strings describing which species to include. By default displays all
            species.
        :type included_species_list: list

        :param return_plotly_figure: Whether or not to return a figure dicctionary of data(graph object traces) and
            layout which may be edited by the user
        :type return_plotly_figure: bool

        :param ddof: Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents
            the number of trajectories. Sample standard deviation uses ddof of 1. Defaults to population standard
            deviation where ddof is 0.
        :type ddof: int

        :param **layout_args: Optional additional arguments to be passed to plotlys layout constructor.
        :type **layout_args: dict
        """

        # Backwards compatibility with xaxis_label argument (which duplicates plotly's xaxis_title argument)
        if layout_args.get('xaxis_title') is not None:
            xaxis_label = layout_args.get('xaxis_title')
            layout_args.pop('xaxis_title')
        if layout_args.get('yaxis_title') is not None:
            yaxis_label = layout_args.get('yaxis_title')
            layout_args.pop('yaxis_title')

        average_trajectory = self.average_ensemble().data[0]
        stddev_trajectory = self.stddev_ensemble(ddof=ddof).data[0]
        from plotly.offline import init_notebook_mode, iplot # pylint: disable=import-outside-toplevel
        import plotly.graph_objs as go # pylint: disable=import-outside-toplevel

        init_notebook_mode(connected=True)

        if not show_title:
            title = 'Mean and Standard Deviation'
        else:
            if title is None:
                title = (self._validate_title(show_title) + " - Mean and Standard Deviation")

        trace_list = []
        for species in average_trajectory:
            if species != 'time':

                if species not in included_species_list and included_species_list:
                    continue

                upper_bound = []
                lower_bound = []
                for i in range(0, len(average_trajectory[species])):
                    upper_bound.append(average_trajectory[species][i] + stddev_trajectory[species][i])
                    lower_bound.append(average_trajectory[species][i] - stddev_trajectory[species][i])

                # Append upper_bound list to trace_list
                trace_list.append(
                    go.Scatter(
                        name=species + ' Upper Bound',
                        x=average_trajectory['time'],
                        y=upper_bound,
                        mode='lines',
                        marker=dict(color="#444"),
                        line=dict(width=1, dash='dot'),
                        legendgroup=str(average_trajectory[species]),
                        showlegend=False
                    )
                )
                trace_list.append(
                    go.Scatter(
                        x=average_trajectory['time'],
                        y=average_trajectory[species],
                        name=species,
                        fillcolor='rgba(68, 68, 68, 0.2)',
                        fill='tonexty',
                        legendgroup=str(average_trajectory[species]),
                    )
                )

                # Append lower_bound list to trace_list
                trace_list.append(
                    go.Scatter(
                        name=species + ' Lower Bound',
                        x=average_trajectory['time'],
                        y= lower_bound,
                        mode='lines',
                        marker=dict(color="#444"),
                        line=dict(width=1, dash='dot'),
                        fillcolor='rgba(68, 68, 68, 0.2)',
                        fill='tonexty',
                        legendgroup=str(average_trajectory[species]),
                        showlegend=False
                    )
                )

        layout = go.Layout(
            showlegend=show_legend,
            title=title,
            xaxis_title=xaxis_label,
            yaxis_title=yaxis_label,
            legend={'traceorder': 'normal'},
            **layout_args
        )
        fig = dict(data=trace_list, layout=layout)

        if not return_plotly_figure:
            iplot(fig)
            return None
        return fig

    def plot_mean_stdev(self, xscale='linear', yscale='linear', xaxis_label="Time", yaxis_label="Value"
                           , title=None, show_title=False, style="default", show_legend=True, included_species_list=[],
                           ddof=0, save_png=False, figsize=(18, 10)):
        """
        Plot a matplotlib graph depicting mean and standard deviation of a results object.

        :param xaxis_label: The label for the x-axis
        :type xaxis_label: str

        :param yaxis_label: The label for the y-axis
        :type yaxis_label: str

        :param title: The title of the graph
        :type title: str

        :param show_title: Default False, if True, title will be displayed on the graph.
        :type show_title: bool

        :param style: Matplotlib style to be displayed on graph.
        :type style: str

        :param show_legend: Default to True, if False, legend will not be shown on graph.
        :type show_legend: bool

        :param included_species_list: A list of strings describing which species to include. By default displays all
            species.
        :type included_species_list: list

        :param ddof: Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents
            the number of trajectories. Sample standard deviation uses ddof of 1. Defaults to population standard
            deviation where ddof is 0.
        :type ddof: int

        :type save_png: bool or str

        :param figsize: The size of the graph. A tuple of the form (width,height). Is (18,10) by default.
        :type figsize: tuple of ints (x,y)
        """

        average_result = self.average_ensemble().data[0]
        stddev_trajectory = self.stddev_ensemble(ddof=ddof).data[0]

        import matplotlib.pyplot as plt # pylint: disable=import-outside-toplevel

        try:
            plt.style.use(style)
        except Exception:
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            msg = f"Invalid matplotlib style. Try using one of the following {plt.style.available}"
            log.warning(msg)
            plt.style.use("default")

        plt.figure(figsize=figsize)

        for species in average_result:
            if species == 'time':
                continue

            if species not in included_species_list and included_species_list:
                continue

            lower_bound = [a-b for a, b in zip(average_result[species], stddev_trajectory[species])]
            upper_bound = [a+b for a, b in zip(average_result[species], stddev_trajectory[species])]

            plt.fill_between(average_result['time'], lower_bound, upper_bound, color='whitesmoke')
            plt.plot(average_result['time'], upper_bound, color='grey', linestyle='dashed')
            plt.plot(average_result['time'], lower_bound, color='grey', linestyle='dashed')
            plt.plot(average_result['time'], average_result[species], label=species)

        if not show_title:
            title = 'Mean and Standard Deviation'
        else:
            if title is None:
                title = (self._validate_title(show_title) + " - Mean and Standard Deviation")

        plt.title(title, fontsize=18)
        plt.xlabel(xaxis_label)
        plt.ylabel(yaxis_label)
        plt.xscale(xscale)
        plt.yscale(yscale)

        plt.plot([0], [11])
        if show_legend:
            plt.legend(loc='best')

        if isinstance(save_png, str):
            plt.savefig(save_png)

        elif save_png:
            plt.savefig(title)

    # for backwards compatability, we need to keep the old name around
    def plotplotly_std_dev_range(self, **kwargs):
        """
        Plot a plotly graph depicting the mean and standard deviation of a results object

        :param xaxis_label: The label for the x-axis
        :type xaxis_label: str

        :param yaxis_label: The label for the y-axis
        :type yaxis_label: str

        :param title: The title of the graph
        :type title: str

        :param show_title: If True, title will be shown on graph.
        :type show_title: bool

        :param show_legend: Default True, if False, legend will not be shown on graph.
        :type show_legend: bool

        :param included_species_list: A list of strings describing which species to include. By default displays all
            species.
        :type included_species_list: list

        :param return_plotly_figure: Whether or not to return a figure dicctionary of data(graph object traces) and
            layout which may be edited by the user
        :type return_plotly_figure: bool

        :param ddof: Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents
            the number of trajectories. Sample standard deviation uses ddof of 1. Defaults to population standard
            deviation where ddof is 0.
        :type ddof: int

        :param **layout_args: Optional additional arguments to be passed to plotlys layout constructor.
        :type **layout_args: dict
        """
        from gillespy2.core import log # pylint: disable=import-outside-toplevel
        log.warning(
            "The plotplotly_std_dev_range function has been deprecated. This function will be removed in a "
            "future release. Please use plotplotly_mean_stdev instead."
        )
        self.plotplotly_mean_stdev(**kwargs)

    # for backwards compatability, we need to keep the old name around
    def plot_std_dev_range(self, **kwargs):
        """
        Plot a matplotlib graph depicting mean and standard deviation of a results object.

        :param xaxis_label: The label for the x-axis
        :type xaxis_label: str

        :param yaxis_label: The label for the y-axis
        :type yaxis_label: str

        :param title: The title of the graph
        :type title: str

        :param show_title: Default False, if True, title will be displayed on the graph.
        :type show_title: bool

        :param style: Matplotlib style to be displayed on graph.
        :type style: str

        :param show_legend: Default to True, if False, legend will not be shown on graph.
        :type show_legend: bool

        :param included_species_list: A list of strings describing which species to include. By default displays all
            species.
        :type included_species_list: list

        :param ddof: Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents
            the number of trajectories. Sample standard deviation uses ddof of 1. Defaults to population standard
            deviation where ddof is 0.
        :type ddof: int

        :type save_png: bool or str

        :param figsize: The size of the graph. A tuple of the form (width,height). Is (18,10) by default.
        :type figsize: tuple of ints (x,y)
        """
        from gillespy2.core import log # pylint: disable=import-outside-toplevel
        log.warning(
            "The plot_std_dev_range function has been deprecated. This function will be removed in a "
            "future release. Please use plot_mean_stdev instead."
        )
        self.plot_mean_stdev(**kwargs)
