import threading
from gillespy2.core import log


class RepeatTimer(threading.Timer):
    """
    Threading timer which repeatedly calls the given function instead of simply ending
    """
    pause = False

    def run(self):
        from IPython.display import clear_output
        type = str.join('', [*self.args[3]])
        self.args = self.args[:3]
        while not self.finished.wait(self.interval):
            self.function(*self.args, **self.kwargs)

        if not self.pause:
            if type == 'progress':
                clear_output()
                print('progress = 100 %')



def display_types():
    return ["graph", "text", "progress"]


def valid_graph_params(live_output_options):
    if live_output_options['type'] not in ['progress', 'graph', 'text']:
        from gillespy2.core.gillespyError import SimulationError
        raise SimulationError("Invalid input to 'live_output', please check spelling and ensure input is"
                              " lower case.")
    if 'interval' not in live_output_options:
        live_output_options['interval'] = 1
    elif live_output_options['interval'] < 0:

        log.warning("In LiveGraphing live_output_options, got \"interval\" = \"{0}\". setting interval = 1"
                    .format(live_output_options['interval']))
        live_output_options['interval'] = 1

    if live_output_options['type'] == "graph" and live_output_options['interval'] < 1:
        log.warning("In LiveGraphing live_output_options, got \"interval\" = \"{0}\". Consider using an interval >= 1 "
                    "when displaying graphs".format(live_output_options['interval']))

    if 'clear_output' not in live_output_options:

        if live_output_options['type'] == "graph" or live_output_options['type'] == "progress":
            live_output_options['clear_output'] = True
        else:
            live_output_options['clear_output'] = False


class LiveDisplayer():
    """
    holds information required for displaying information when live_output = True
    """

    def __init__(self, model=None, timeline=None, number_of_trajectories=1, live_output_options={}, resume=False):

        self.display_type = live_output_options['type']
        self.display_interval = live_output_options['interval']
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
        return "Trajectory (" + str(self.current_trajectory) + "/" + str(self.number_of_trajectories) + ")"

    def increment_trajectory(self, trajectory_num):
        self.current_trajectory = trajectory_num + 1
        self.header_printed = False

    def print_text_header(self):

        self.header_printed = True
        if self.number_of_trajectories > 1:
            print(self.trajectory_header())

        print("Time      |", end="")
        for species in self.model.listOfSpecies:
            print(species[:10].ljust(10), end="|")
        print("")

    '''
    curr_state and curr_time should be list of len 1 to get reference
    '''
    def display(self, curr_state, curr_time, trajectory_base):
        from IPython.display import clear_output
        from math import floor
        curr_time = curr_time[0]
        curr_state = curr_state[0]

        # necessary for __f function in hybrid solver
        if 't' in curr_state:
            if curr_state['t'] > curr_time:
                curr_time = curr_state['t']
        elif 'time' in curr_state:
            if curr_state['time'] > curr_time:
                curr_time = curr_state['time']

        if self.clear_output:
                clear_output(wait=True)

        if self.display_type == "text":

            if not self.header_printed:
                self.print_text_header()

            print(str(round(curr_time, 2))[:10].ljust(10), end="|")

            for i in range(self.number_species):
                print(str(curr_state[self.species[i]])[:10].ljust(10), end="|")
            print("")

        elif self.display_type == "progress":

            if self.number_of_trajectories > 1:
                print(self.trajectory_header())
            if self.resume is True:
                print("progress =", round(((curr_time-self.x_shift)/self.timeline_len)*100, 2), "%\n"
                      )
            else:
                print("progress =", round((curr_time / (self.timeline_len + self.x_shift)) * 100, 2), "%\n")

        elif self.display_type == "graph":

            import matplotlib.pyplot as plt
            from gillespy2.core.results import common_rgb_values

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
