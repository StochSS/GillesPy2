from gillespy2.core import log

import threading
class RepeatTimer(threading.Timer):
    def run(self):
        while not self.finished.wait(self.interval):
            self.function(*self.args, **self.kwargs)

def display_types():
    return ["graph","text","progress"]

def valid_graph_params(display_type,display_interval):

    if display_interval < 0:
        log.warning("Got display_interval = \"{0}\". Must use display_interval > 0".format(display_interval))
        return False

    elif display_interval > 0:

        if display_type in display_types():

            if display_type == "graph" and display_interval < 1:
                log.warning("Got display_interval = \"{0}\". Consider using a longer interval when displaying graphs"
                            .format(display_interval))

            return True

        #display_type will default to progress if "None"
        elif display_type == None:
            return True


        else:
            log.warning(
                'Got display_type = \"{0}\". Display_type should be \"graph\", \"text\", or \"progress\"'.format(
                    display_type))

            return False

class LiveDisplayer():

    def __init__(self,display_type = None,display_interval = 0,model = None,timeline_len=None,number_of_trajectories=1):

        self.display_type = display_type
        self.display_interval = display_interval
        self.model = model
        self.timeline_len = timeline_len
        self.number_of_trajectories = number_of_trajectories

        species_mappings = model._listOfSpecies
        self.species = list(species_mappings.keys())

        self.number_species = len(self.species)
        self.current_trajectory = 1
        self.header_printed = False

        if display_type is None:
                self.display_type = "progress"
                log.warning('Unspecified display_type. Displaying progress.')

    def trajectory_header(self):
        return "Trajectory ("+ str(self.current_trajectory)+ "/"+ str(self.number_of_trajectories)+ ")"

    def increment_trajectory(self,trajectory_num):
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

        try:
            if self.display_type == "text":

                if not self.header_printed:
                    self.print_text_header()

                print(str(round(curr_time, 2))[:10].ljust(10), end="|")

                for i in range(self.number_species):
                    print(str(curr_state[self.species[i]])[:10].ljust(10), end="|")
                print("")

            elif self.display_type == "progress":

                clear_output(wait=True)
                if self.number_of_trajectories > 1:
                    print(self.trajectory_header())

                print("progress =", round((curr_time / self.timeline_len) * 100, 2), "%\n")

            elif self.display_type == "graph":

                import matplotlib.pyplot as plt
                from gillespy2.core.results import common_rgb_values

                entry_count = floor(curr_time)

                clear_output(wait=True)

                plt.figure(figsize=(18, 10))
                plt.xlim(right=self.timeline_len)
                plt.title(self.trajectory_header())

                for i in range(self.number_species):
                    line_color = common_rgb_values()[(i) % len(common_rgb_values())]

                    plt.plot(trajectory_base[0][:, 0][:entry_count].tolist(),
                             trajectory_base[0][:, i + 1][:entry_count].tolist(), color=line_color,
                             label=self.species[i])

                    plt.plot([entry_count - 1, curr_time], [trajectory_base[0][:, i + 1][entry_count - 1],
                                                            curr_state[self.species[i]]], linewidth=3,
                             color=line_color)

                plt.legend(loc='upper right')
                plt.show()
        except:
            log.warning("exception in liveGraphing.display. Variables may not have initialized properly before display was called.")