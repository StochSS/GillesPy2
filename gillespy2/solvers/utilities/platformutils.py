"""
Collection of utility functions for performing platform-specific simulations.
Necessary to accomodate for things like signals and processes which vary between
Windows and POSIX platforms.

Written By: Joshua Cooper
April 18, 2021
"""

import subprocess
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
