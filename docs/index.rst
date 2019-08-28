Documentation for GillesPy2 |release|
#####################################

GillesPy2 is an open-source Python package for stochastic simulation of biochemical systems.  It offers an object-oriented approach for creating mathematical models of biological systems, as well as a variety of methods for performing time simulation of those models.  The methods include the `Gillespie direct method (SSA) <https://en.wikipedia.org/wiki/Gillespie_algorithm>`_, several variant stochastic simulation methods including `tau leaping <https://en.wikipedia.org/wiki/Tau-leaping>`_, and numerical integration of ODEs.  The solvers support a variety of user environments, with optimized code for C++, `Cython <https://cython.org>`_, and `NumPy <https://numpy.org>`_.  Models can also be read from files in `SBML <https://en.wikipedia.org/wiki/SBML>`_ format.


Getting a copy of GillesPy2
***************************

The latest version of GillesPy2 can be found on `PyPI <https://pypi.org/project/gillespy2>`_.  The source code is available on `GitHub <https://github.com/GillesPy2/GillesPy2>`_.  GillesPy2 is licensed under the GNU General Public License version 3.

.. raw:: html
   
   <p style="width: 80%; margin: auto; padding: 0.5em; border: 1px solid purple"><a href="https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform"><b>Please register as a user!</b></a>  GillesPy2's development is funded by NIH grant 2R01EB014877, and to continue support, we need to prove GillesPy2 has users. <a href="https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform">Please fill out our short registration form!</a></p>



Getting help
************

If you find any problem with GillesPy2 or this documentation, please report it using `the GitHub issue tracker <https://github.com/GillesPy2/GillesPy2/issues>`_ for the project.  You can also contact the main author, Dr. `Brian Drawert <http://www.cs.unca.edu/~drawert>`_, directly with questions and suggestions.


Documentation
*************

.. toctree::
   :maxdepth: 1
   :caption: Getting started
   :name: sec-getting-started

   getting_started/installation/installation.rst
   getting_started/basic_usage/basic_usage.rst


.. toctree::
   :maxdepth: 1
   :caption: Tutorials
   :name: sec-tutorials

   tutorials/tut_toggle_switch/tut_toggle_switch.rst
   tutorials/tut_michaelis_menten/tut_michaelis_menten.rst
   tutorials/tut_sbml/tut_sbml.rst


.. toctree::
   :maxdepth: 3
   :caption: API reference

   classes/gillespy2


Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
