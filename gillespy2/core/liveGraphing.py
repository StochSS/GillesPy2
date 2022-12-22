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

import sys
import json
import threading

from math import floor

from gillespy2.core.results import common_rgb_values
from gillespy2.core.gillespyError import SimulationError
from gillespy2.core import log

class CRepeatTimer(threading.Timer):
    """
    Threading timer which repeatedly calls the given function instead of simply ending.

    Used for C solver live graphing support
    """
    pause = False

    def run(self):
        _ = str.join('', [*self.args[1]])
        while not self.finished.wait(self.interval):
            args = self.args[0].get()
            self.function(*args, **self.kwargs)

        if not self.pause:
            args = self.args[0].get()
            self.kwargs['finished'] = True
            self.function(*args, **self.kwargs)


class RepeatTimer(threading.Timer):
    """
    Threading timer which repeatedly calls the given function instead of simply ending
    """
    pause = False

    def run(self):
        _ = str.join('', [*self.args[3]])
        self.args = self.args[:3]
        while not self.finished.wait(self.interval):
            self.function(*self.args, **self.kwargs)

        if not self.pause:
            self.kwargs['finished'] = True
            self.function(*self.args, **self.kwargs)



def display_types():
    '''
    Get the list of supported display types.

    :returns: Supported display types.
    :rtype: list
    '''
    return ["graph", "text", "progress"]


def valid_graph_params(live_output_options):
    '''
    Validated the live output options.

    :param live_output_options: Options to be validated.
    :type live_output_options: dict

    :raises SimulationError: If the display type is invalid.
    '''
    if live_output_options['type'] not in ['progress', 'graph', 'text']:
        raise SimulationError("Invalid input to 'live_output', please check spelling and ensure input is"
                              " lower case.")
    if 'interval' not in live_output_options:
        live_output_options['interval'] = 1
    elif live_output_options['interval'] < 0:
        message = f"In LiveGraphing live_output_options, got 'interval' = '{live_output_options['interval']}'."
        message += " setting interval = 1"
        log.warning(message)
        live_output_options['interval'] = 1

    if live_output_options['type'] == "graph" and live_output_options['interval'] < 1:
        message = f"In LiveGraphing live_output_options, got 'interval' = '{live_output_options['interval']}'."
        message += "Consider using an interval >= 1 when displaying graphs."
        log.warning(message)

    if 'clear_output' not in live_output_options:
        if live_output_options['type'] == "graph" or live_output_options['type'] == "progress":
            live_output_options['clear_output'] = True
        else:
            live_output_options['clear_output'] = False

    if 'file_path' not in live_output_options:
        live_output_options['file_path'] = None
    elif live_output_options['type'] == "graph" and live_output_options['file_path'] is not None:
        live_output_options['type'] = "figure"


class LiveDisplayer():
    """
    holds information required for displaying information when live_output = True
    """

    def __init__(self, model=None, timeline=None, number_of_trajectories=1, live_output_options={}, resume=False):

        self.display_type = live_output_options['type']
        self.display_interval = live_output_options['interval']
        self.file_path = live_output_options['file_path']
        self.model = model
        self.resume = resume
        self.timeline = timeline
        self.timeline_len = timeline.size
        self.x_shift = int(timeline[0])
        self.number_of_trajectories = number_of_trajectories
        self.clear_output = live_output_options['clear_output']
        species_mappings = model._listOfSpecies
        self.species = list(species_mappings.keys())
        self.number_species = len(self.species)
        self.current_trajectory = 1
        self.header_printed = False

    def trajectory_header(self):
        '''
        Create the trajectory header for the output.
        '''
        return "Trajectory (" + str(self.current_trajectory) + "/" + str(self.number_of_trajectories) + ")"

    def increment_trajectory(self, trajectory_num):
        '''
        Increment the trejectory counter.
        '''
        self.current_trajectory = trajectory_num + 1
        self.header_printed = False

    def print_text_header(self, file_obj):
        '''
        Print the header for text display type.

        :param file_obj: File object to write text output.
        :type file_obj: file object
        '''
        self.header_printed = True
        if self.number_of_trajectories > 1:
            print(self.trajectory_header(), file=file_obj)

        print("Time      |", end="", file=file_obj)
        for species in self.model.listOfSpecies:
            print(species[:10].ljust(10), end="|", file=file_obj)
        print("", file=file_obj)

    def display(self, curr_state, curr_time, trajectory_base, finished=False):
        '''
        Display the output for the live grapher.

        :param curr_state: Current state of the simulation. Should be a list of len 1 to get reference.
        :type curr_state: list

        :param curr_time: Current time of the simulation. Should be a list of len 1 to get reference.
        :type curr_time: list

        :param trajectory_base: Current results of the simulation.
        :type trajectory_base: list

        :param finished: Indicates whether or not the simulation has finished.
        :type finished: bool
        '''
        from IPython.display import clear_output # pylint: disable=import-outside-toplevel

        curr_time = curr_time[0]
        curr_state = curr_state[0]

        # necessary for __f function in hybrid solver
        if 't' in curr_state:
            if curr_state['t'] > curr_time:
                curr_time = curr_state['t']
        elif 'time' in curr_state:
            if curr_state['time'] > curr_time:
                curr_time = curr_state['time']

        if self.file_path is None or self.display_type == "graph":
            if self.clear_output:
                clear_output(wait=True)
            file_obj = sys.stdout
        else:
            mode = "w" if self.clear_output else "a"
            file_obj = open(self.file_path, mode, encoding="utf-8")

        if self.display_type == "text":

            if not self.header_printed:
                self.print_text_header(file_obj)

            print(str(round(curr_time, 2))[:10].ljust(10), end="|", file=file_obj)

            for i in range(self.number_species):
                print(str(curr_state[self.species[i]])[:10].ljust(10), end="|", file=file_obj)
            print("", file=file_obj)

        elif finished and self.display_type == "progress":
            print("progress = 100 %", file=file_obj)

        elif self.display_type == "progress":
            if self.number_of_trajectories > 1:
                print(self.trajectory_header(), file=file_obj)
            if self.resume is True:
                print(f"progress = {round(((curr_time-self.x_shift)/self.timeline_len)*100, 2)} %\n", file=file_obj)
            else:
                print(
                    f"progress = {round((curr_time / (self.timeline_len + self.x_shift)) * 100, 2) }%\n", file=file_obj
                )

        elif self.display_type == "graph":
            if finished:
                return

            import matplotlib.pyplot as plt # pylint: disable=import-outside-toplevel

            entry_count = floor(curr_time) - self.x_shift

            plt.figure(figsize=(18, 10))
            plt.xlim(right=self.timeline[-1])
            plt.xlim(left=self.timeline[0])
            plt.title(self.trajectory_header())
            for i in range(self.number_species):
                line_color = common_rgb_values()[(i) % len(common_rgb_values())]
                plt.plot(trajectory_base[0][:, 0][:entry_count].tolist(),
                         trajectory_base[0][:, i + 1][:entry_count].tolist(), color=line_color,
                         label=self.species[i])

                plt.plot([entry_count - 1, curr_time - self.timeline[0]], [trajectory_base[0][:, i + 1][entry_count - 1]
                    , curr_state[self.species[i]]], linewidth=3,
                         color=line_color)
            plt.legend(loc='upper right')
            plt.show()

        elif self.display_type == "figure":
            import plotly # pylint: disable=import-outside-toplevel
            import plotly.graph_objs as go # pylint: disable=import-outside-toplevel

            entry_count = floor(curr_time) - self.x_shift

            trace_list = []
            for i, species in enumerate(self.species):
                line_dict = {"color": common_rgb_values()[(i) % len(common_rgb_values())]}
                trace_list.append(
                    go.Scatter(
                        x=trajectory_base[0][:, 0][:entry_count].tolist(),
                        y=trajectory_base[0][:, i + 1][:entry_count].tolist(),
                        mode="lines", name=species, line=line_dict, legendgroup=species
                    )
                )

            layout = go.Layout(
                showlegend=True, title=self.trajectory_header(),
                xaxis={"range": [self.timeline[0], self.timeline[-1]]}
            )

            fig = dict(data=trace_list, layout=layout)
            json.dump(fig, file_obj, cls=plotly.utils.PlotlyJSONEncoder)

        if self.file_path is not None and self.display_type != "graph":
            file_obj.close()
