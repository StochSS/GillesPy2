#!/usr/bin/env python3
# =============================================================================
# @file    setup.py
# @brief   GillesPy2 setup file
# @license Please see the file named LICENSE in the project directory
# @website https://github.com/GillesPy2/GillesPy2
#
# Note: how to do a PyPI release
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run the following commands:
#
#   python3 setup.py sdist bdist_wheel
#   twine upload dist/*
#
# =============================================================================

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.easy_install import easy_install
import os
from os import path


# Read the contents of auxiliary files.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SETUP_DIR = path.dirname(os.path.abspath(__file__))

with open(path.join(SETUP_DIR, 'requirements.txt')) as f:
    reqs = f.read().rstrip().splitlines()

with open(path.join(SETUP_DIR, 'README.md'), 'r', errors = 'ignore') as f:
    readme = f.read()

# The following reads the variables without doing an "import handprint",
# because the latter will cause the python execution environment to fail if
# any dependencies are not already installed -- negating most of the reason
# we're using setup() in the first place.  This code avoids eval, for security.

version = {}
with open(path.join(SETUP_DIR, 'gillespy2/__version__.py')) as f:
    text = f.read().rstrip().splitlines()
    vars = [line for line in text if line.startswith('__') and '=' in line]
    for v in vars:
        setting = v.split('=')
        version[setting[0].strip()] = setting[1].strip().replace("'", '')


# Decorate some command classes used below.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def stoch_path(command_subclass):
    """
    A decorator for classes subclassing one of the setuptools commands.
    It modifies the run() method.
    """
    orig_run = command_subclass.run

    def modified_run(self):
        success = False
        orig_run(self)

    command_subclass.run = modified_run
    return command_subclass


@stoch_path
class develop_new(develop):
    pass


@stoch_path
class install_new(install):
    pass


@stoch_path
class bdist_egg_new(bdist_egg):
    pass


@stoch_path
class easy_install_new(easy_install):
    pass


# Finally, define our namesake.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setup(name             = version['__title__'].lower(),
      version          = version['__version__'],
      author           = version['__author__'],
      author_email     = version['__email__'],
      maintainer       = version['__author__'],
      maintainer_email = version['__email__'],
      license          = version['__license__'],
      url              = version['__url__'],
      download_url     = version['__download_url__'],
      description      = version['__description__'],
      long_description = readme,
      long_description_content_type = "text/markdown",
      keywords         = "gillespy2, Gillespie algorithm, biochemical simulation",
      project_urls     = {
          "Tracker": "https://github.com/GillesPy2/GillesPy2/issues",
          "Source" : "https://github.com/GillesPy2/GillesPy2",
      },

      packages         = find_packages('.'),
      include_package_data = True,
      install_requires = reqs,
      cmdclass         = {'bdist_egg': bdist_egg_new,
                          'install': install_new,
                          'develop': develop_new,
                          'easy_install': easy_install_new},

      classifiers      = [
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Intended Audience :: Science/Research'
      ],
)
