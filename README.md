GillesPy2
=========

GillesPy2 is a Python package for stochastic simulation of biochemical systems.  It offers an object-oriented approach for creating mathematical models of biological systems, as well as a variety of methods for performing time simulation of those models.  The methods include the [Gillespie direct method (SSA)](https://en.wikipedia.org/wiki/Gillespie_algorithm), several variant stochastic simulation methods including [tau-leaping](https://en.wikipedia.org/wiki/Tau-leaping), and numerical integration of ODEs.  The solvers support a variety of runtime environments, with optimized code for C++, [Cython](https://cython.org), and [NumPy](https://numpy.org).  GillesPy2 also supports [SBML](https://en.wikipedia.org/wiki/SBML).

<table><tr><td><b>
<img width="20%" align="right" src=".graphics/stochss-logo.svg">
GillesPy2 is part of the <a href="http://www.stochss.org">StochSS</a> project. We rely on grants to sustain development.<br><a href="https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform">PLEASE REGISTER AS A USER</a>, so that we can demonstrate GillesPy2's popularity!</b>
</td></tr></table>

![PyPI](https://img.shields.io/pypi/v/gillespy2.svg)
![PyPI - License](https://img.shields.io/pypi/l/gillespy2.svg?color=green)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/gillespy2.svg)


## Installation
You can install gillespy2 through the following methods.

### Linux

```
python3 -m pip install --upgrade pip
python3 -m pip install gillespy2
```
### Mac
(Pip should work, same as Linux)

### Windows
(Pip should work, same as Linux)


## Getting Started

GillesPy2 performs stochastic biochemical simulations utilizing simplified object-oriented abstractions.  The basic steps to use GillesPy2 are:

1. create (or SBML import) a GillesPy2.Model containing species, parameters, and reactions.
2. call the model's '.run()' method.

The run method can be customized using keyword arguments to select different solvers, random seed, data return type and more.
For more detailed examples on how to use GillesPy2, please see the [Getting Started](https://github.com/GillesPy2/GillesPy2/tree/master/examples/Getting-Started.ipynb) notebook contained in the [examples](https://github.com/GillesPy2/GillesPy2/tree/master/examples) subdirectory.

## Built With

* [Numpy](http://www.numpy.org/)
* [Scipy](https://www.scipy.org/)
* [MatplotLib](https://matplotlib.org/index.html)

## Contributing

If you have any problems, or want to request a feature, please submit an issue to this repository.  If you want to contribute to GillesPy2, please follow the guidelines set forth in [CONTRIBUTING.md](https://github.com/GillesPy2/GillesPy2/tree/master/CONTRIBUTING.md).  If you have any questions, contact Brian Drawert.

## Authors

* **Dr. Brian Drawert** 
* **Dr. Kevin Sanft**
* **Ghilman Brock**  
* **Eliot Dixon**  
* **Dalton Nickerson**  
* **Sean Matthew**
* **George Hall** 
* **W.R. Jackson** 
* **Samuel Hodges**
* **Emma Weisgerber**

## License

GillesPy2 is licenced under GPLv3, see [LICENCE] for details.


## Acknowledgments
This work has been funded by National Institutes of Health (NIH) NIBIB Award No. 2R01EB014877-04A1

## Build Status

| Master Branch |  Develop Branch | Coverage |
|----------------|---|---|
| [![Build Status](https://travis-ci.org/GillesPy2/GillesPy2.svg?branch=master)](https://travis-ci.org/GillesPy2/GillesPy2) | [![Build Status](https://travis-ci.org/GillesPy2/GillesPy2.svg?branch=develop)](https://travis-ci.org/GillesPy2/GillesPy2) | ![Coverage](coverage.svg) |
