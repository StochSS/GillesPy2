from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.easy_install import easy_install
import os

SETUP_DIR = os.path.dirname(os.path.abspath(__file__))


def stoch_path(command_subclass):
    """
    A decorator for classes subclassing one of the setuptools commands.
    It modifies the run() method.
    """
    orig_run = command_subclass.run
    
    def modified_run(self):
      
        success=False
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


setup(name = "gillespy2",
      version = "1.1",
      packages = find_packages('.'),
      include_package_data=True,
      description = 'Python interface for Gillespie style biochemical simulations',
      
      install_requires = ["numpy",
                          "matplotlib",
                          "scipy"],
      
      author = "Brian Drawert, Kevin Sanft, Ghilman Brock, Eliot Dixon, Dalton Nickerson",
      author_email = ["bdrawert@unca.edu"],
      license = "GPL",
      keywords = "gillespy2, gillespie algorithm, biochemical simulation",

      url = "http://www.github.com/briandrawert/GillesPy2", # we don't really yet have one

      download_url = "https://github.com/briandrawert/GillesPy2/tarball/master/",
      
      cmdclass = {'bdist_egg':bdist_egg_new,
                  'install':install_new,
                  'develop':develop_new,
                  'easy_install':easy_install_new}
      )
