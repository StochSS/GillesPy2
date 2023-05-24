//
// Created by josh on 2/28/2023.
//

#include "model.h"
#include "simulation.h"
#include "extension.h"
#include "ssa_cpp_solver/SSASimulation.h"

#include <vector>

static PyMethodDef cgillespy_module_methods[] = {
        {
            "solve_ssa",
            (PyCFunction) Gillespy::solve_ssa,
            METH_KEYWORDS | METH_VARARGS,
            "Run the GillesPy2 model using a compiled SSA solver implementation"
        },
        {
                "solve_ode",
                (PyCFunction) Gillespy::solve_ode,
                METH_KEYWORDS | METH_VARARGS,
                "Run the GillesPy2 model using a compiled ODE solver implementation"
        },
        {
                "solve_tau_leaping",
                (PyCFunction) Gillespy::solve_tau_leaping,
                METH_KEYWORDS | METH_VARARGS,
                "Run the GillesPy2 model using a compiled Tau-Leaping solver implementation"
        },
        {
                "solve_tau_hybrid",
                (PyCFunction) Gillespy::solve_tau_hybrid,
                METH_KEYWORDS | METH_VARARGS,
                "Run the GillesPy2 model using a compiled Tau-Hybrid solver implementation"
        },
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cgillespy_module_def = {
        PyModuleDef_HEAD_INIT,
        "libcgillespy",
        NULL,
        -1,
        cgillespy_module_methods
};

static PyObject *Gillespy::solve_ode(PyObject *self, PyObject *args, PyObject *kwargs)
{
    return PyUnicode_FromString("not implemented");
}

static PyObject *Gillespy::solve_ssa(PyObject *self, PyObject *args, PyObject *kwargs)
{
    std::vector<char*> argv;
    static char *keywords[] = {
            (char*)"argv",
            NULL
    };
    PyObject *py_argv;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", keywords, &py_argv))
    {
        std::cerr << "bad arguments" << std::endl;
        return nullptr;
    }
    int list_size = PyList_Size(py_argv);

    for (int list_index = 0; list_index < list_size; ++list_index)
    {
        PyObject *py_arg = PyList_GetItem(py_argv, list_index);
        char *arg = PyBytes_AsString(PyUnicode_AsUTF8String(py_arg));
        argv.emplace_back(arg);
    }

    long return_code = Gillespy::run_ssa_simulation(list_size, argv.data());
    return PyLong_FromLong(return_code);
}

static PyObject *Gillespy::solve_tau_leaping(PyObject *self, PyObject *args, PyObject *kwargs)
{
    return PyUnicode_FromString("not implemented");
}

static PyObject *Gillespy::solve_tau_hybrid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    return PyUnicode_FromString("not implemented");
}

PyMODINIT_FUNC
PyInit_libcgillespy(void)
{
    return PyModule_Create(&cgillespy_module_def);
}
