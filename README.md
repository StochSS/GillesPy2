# GillesPy2

  GillesPy2 is a python library designed to create stochastic simulations for biochemical systems.  This package provides an object-oriented approach to creating mathematical models based on these real-world systems and simulates reaction events over time, selecting from a variety of algorithms including ODE solutions, the Gillespie Direct algorithm (SSA), and several variants of the SSA.  This library contains multiple versions of solvers to support a variety of user environments with optimized code for C++, Cython, and NumPy.  
  
**GillesPy2 is part of the StochSS project [http://www.stochss.org/], and we are relying on continued funding for sustained development. Please consider registering to show your support. Register [here](https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform):**  
  
https://docs.google.com/forms/d/12tAH4f8CJ-3F-lK44Q9uQHFio_mGoK0oY829q5lD7i4/viewform

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
For more detailed examples on how to use GillesPy2, please see the [Getting Started](/examples/Getting-Started.ipynb) notebook contained in the [examples](/examples) subdirectory.

## Built With

* [Numpy](http://www.numpy.org/)
* [Scipy](https://www.scipy.org/)
* [MatplotLib](https://matplotlib.org/index.html)

## Contributing

If you have any problems, or want to request a feature, please submit an issue to this repository.  If you want to contribute to GillesPy2, please follow the guidelines set forth in [CONTRIBUTING.md](CONTRIBUTING.md).  If you have any questions, contact Brian Drawert.

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
