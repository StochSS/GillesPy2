from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.easy_install import easy_install
import os

##########################################################
### HOW TO DO A PYPI RELEASE ###
# python3 setup.py sdist bdist_wheel
# twine upload dist/* 
##########################################################

SETUP_DIR = os.path.dirname(os.path.abspath(__file__))


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


# update all install classes with our new class
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


with open('README.md', 'r') as fh:
    full_description = fh.read()

setup(name="gillespy2",
      version="1.1.2",
      packages=find_packages('.'),
      include_package_data=True,
      description='Python interface for Gillespie style biochemical simulations',
      long_description=full_description,
      long_description_content_type="text/markdown",

      install_requires=["numpy",
                        "matplotlib",
                        "scipy"],

      author="Brian Drawert, Kevin Sanft, Sean Matthew, George Hall, Dalton Nickerson, Samuel Hodges, Emma Weisgerber, Eliot Dixon, Ghilman Brock, W.R. Jackson",
      author_email="bdrawert@unca.edu",
      license="GPL",
      keywords="gillespy2, gillespie algorithm, biochemical simulation",

      url="https://gillespy2.github.io/GillesPy2/",

      download_url="https://github.com/GillesPy2/GillesPy2/tarball/master/",

      cmdclass={'bdist_egg': bdist_egg_new,
                'install': install_new,
                'develop': develop_new,
                'easy_install': easy_install_new},
      classifiers=[
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
