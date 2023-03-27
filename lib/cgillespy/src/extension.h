#pragma once

#include <Python.h>

namespace Gillespy
{
    static PyObject *solve_ssa(PyObject *self, PyObject *args, PyObject *kwargs);
    static PyObject *solve_ode(PyObject *self, PyObject *args, PyObject *kwargs);
    static PyObject *solve_tau_leaping(PyObject *self, PyObject *args, PyObject *kwargs);
    static PyObject *solve_tau_hybrid(PyObject *self, PyObject *args, PyObject *kwargs);
}
