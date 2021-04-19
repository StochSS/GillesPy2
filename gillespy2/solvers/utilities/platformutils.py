"""
Collection of utility functions for performing platform-specific simulations.
Necessary to accomodate for things like signals and processes which vary between
Windows and POSIX platforms.

Written By: Joshua Cooper
April 18, 2021
"""

import subprocess
import signal
import os

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
