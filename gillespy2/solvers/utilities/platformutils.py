"""
Collection of utility functions for performing platform-specific simulations.
Necessary to accomodate for things like signals and processes which vary between
Windows and POSIX platforms.

Written By: Joshua Cooper
April 18, 2021
"""

import subprocess
import threading
import time
import signal
import os

class SimulationReader():
    """
    """
    def __init__(self, simulation: subprocess.Popen, timeout=0):
        self.simulation = simulation
        self.buffer = []
        self.timed_out = False

        if timeout > 0:
            self.timeout_thread = threading.Timer(timeout, self.__timeout_kill)
        else:
            self.timeout_thread = None
        self.reader_thread  = threading.Thread(name="SimulationReaderThread",
                                               target=self.__reader_thread)

    def __timeout_kill(self):
        """
        Handler to kill the simulation on timeout.
        """
        self.timed_out = True
        sub_kill(self.simulation)

    def __reader_thread(self):
        """
        Handler for reading data from subprocess in background thread.
        """
        def read_next():
            # Reads the next block from the simulation output.
            # Returns the length of the string read.
            line = self.simulation.stdout.read().decode("utf-8")
            ln = len(line)
            if ln > 0:
                self.buffer.append(line)
            return ln
        # Read output 1 block at a time, until the program is finished.
        page_size = read_next()
        while page_size > 0 and self.simulation.poll() is None:
            page_size = read_next()

    def start(self):
        """
        Activates the reader by starting the processing/timeout threads.
        Must be called before read().
        """
        if self.timeout_thread is not None:
            self.timeout_thread.start()
        self.reader_thread.start()

    def read(self) -> (str, int):
        """
        Blocks and waits for the output of the simulation to finish.
        Returns the output as a string, and the return code of the simulation.
        """
        return_code = self.simulation.poll()

        while return_code is None:
            return_code = self.simulation.poll()
            time.sleep(0.1)

        self.reader_thread.join()
        if self.timeout_thread is not None:
            self.timeout_thread.cancel()
            self.timeout_thread.join()

        return "".join(self.buffer), return_code


def open_simulation(args, **kwargs):
    """
    A thin wrapper around the Popen class which adds platform-specific args.
    Pass arguments just like you would Popen.
    """
    # Append universally-needed flags
    kwargs["start_new_session"] = True

    # Windows-specific flags
    if os.name == "nt":
        kwargs["creationflags"] |= subprocess.CREATE_NEW_PROCESS_GROUP

    return subprocess.Popen(args, **kwargs)

def sub_kill(simulation: subprocess.Popen):
    """
    Platform-independent wrapper around `os.kill` and similar functions.
    Sends an interrupt signal/event to the given simulation.
    """
    # Windows-specific interrupt
    if os.name == "nt":
        simulation.send_signal(signal.CTRL_BREAK_EVENT)
    # POSIX interrupt
    else:
        os.killpg(simulation.pid, signal.SIGINT)
