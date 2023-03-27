//
// Created by josh on 2/28/2023.
//

#include "model.h"
#include "simulation.h"
#include "extension.h"
#include "ssa_cpp_solver/SSASolver.h"

#include <vector>
#include <string>

static PyObject *model_ex(PyObject *self, PyObject *args, PyObject *kwargs);

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
        {
            "model_ex",
            (PyCFunction) model_ex,
            METH_KEYWORDS | METH_VARARGS,
            "Crash test of model run implementation (do not publish!)"
        },
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cgillespy_module_def = {
        PyModuleDef_HEAD_INIT,
        "cgillespy",
        NULL,
        -1,
        cgillespy_module_methods
};

static PyObject *model_ex(PyObject *self, PyObject *args, PyObject *kwargs)
{
    // assumed structure:
    //    species = dict[str, float] : key is name, float is initial population
    //    reactions = dict[str, list[list[float]]] : key is rx name, 2d list is stoichiometry
    //    parameters = dict[str, float] : key is parameter name, float is value
    static char *keywords[] = {
            // (char*) "model",
            (char*) "species",
            (char*) "reactions",
            (char*) "parameters",
//            (char*) "rate_rules",
//            (char*) "events",
            NULL
    };
    PyObject *py_species;
    PyObject *py_reactions;
    PyObject *py_parameters;
    PyObject *py_rate_rules;
    PyObject *py_events;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO", keywords, &py_species, &py_reactions, &py_parameters))
    {
        return NULL;
    }

    if (!PyDict_Check(py_species) || !PyDict_Check(py_reactions) || !PyDict_Check(py_parameters))
    {
        return NULL;
    }

    std::vector<std::string> species_names;
    std::vector<double> initial_populations;
    std::vector<std::string> reaction_names;
    PyObject *py_dict_key, *py_dict_val;
    Py_ssize_t py_dict_pos;

    py_dict_pos = 0;
    while (PyDict_Next(py_species, &py_dict_pos, &py_dict_key, &py_dict_val))
    {
        const char *species_name = PyBytes_AsString(PyUnicode_AsUTF8String(py_dict_key));
        std::cerr << "species_name = " << species_name << std::endl;
        double initial_value = PyFloat_AsDouble(py_dict_val);
        std::cerr << "intial_value = " << initial_value << std::endl;
        initial_populations.push_back(initial_value);
        species_names.emplace_back(species_name);
    }

    py_dict_pos = 0;
    while (PyDict_Next(py_reactions, &py_dict_pos, &py_dict_key, &py_dict_val))
    {
        const char *reaction_name = PyBytes_AsString(PyUnicode_AsUTF8String(py_dict_key));
        std::cerr << "reaction_name = " << reaction_name << std::endl;
        reaction_names.emplace_back(reaction_name);
    }

    std::cerr << "got here" << std::endl;
    Gillespy::Model<unsigned int> model(species_names, initial_populations, reaction_names);
    model.number_reactions = PyDict_Size(py_reactions);
    py_dict_pos = 0;
    while (PyDict_Next(py_reactions, &py_dict_pos, &py_dict_key, &py_dict_val))
    {
        PyObject *reactants = PyList_GetItem(py_dict_val, 0);
        PyObject *products = PyList_GetItem(py_dict_val, 1);
        std::cerr << "reaction #" << py_dict_pos << std::endl;

        Py_ssize_t list_pos, list_size;
        for (list_pos = 0, list_size = PyList_Size(products); list_pos < list_size; ++list_pos)
        {
            std::cerr << "product (" << list_pos << "/" << list_size << ") = " << std::endl;
            PyObject *py_product = PyList_GetItem(products, list_pos);
            int product = PyLong_AsLong(py_product);
            std::cerr << "\t" << product << std::endl;
            model.reactions[py_dict_pos - 1].products_change[list_pos] = product;
            model.reactions[py_dict_pos - 1].species_change[list_pos]  = product;
        }
        for (list_pos = 0, list_size = PyList_Size(reactants); list_pos < list_size; ++list_pos)
        {
            std::cerr << "reactant (" << list_pos << "/" << list_size << ")" << std::endl;
            PyObject *py_reactant = PyList_GetItem(reactants, list_pos);
            int reactant = PyLong_AsLong(py_reactant);
            std::cerr << "\t" << reactant << std::endl;
            model.reactions[py_dict_pos - 1].reactants_change[list_pos] = reactant;
            model.reactions[py_dict_pos - 1].species_change[list_pos]  -= reactant;
        }
    }
    model.update_affected_reactions();

    Gillespy::Simulation<unsigned int> simulation;
    simulation.output_interval = 1;
    simulation.model = &model;
    simulation.end_time = 20;
    simulation.random_seed = time(0);
    simulation.number_timesteps = 21;
    simulation.number_trajectories = 1;
    Gillespy::init_simulation(&model, simulation);
    Gillespy::ssa_direct(&simulation);
    simulation.output_buffer_final(std::cout);

    Py_RETURN_NONE;
}

static PyObject *Gillespy::solve_ode(PyObject *self, PyObject *args, PyObject *kwargs)
{
    return PyUnicode_FromString("not implemented");
}

static PyObject *Gillespy::solve_ssa(PyObject *self, PyObject *args, PyObject *kwargs)
{
    return PyUnicode_FromString("not implemented");
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
PyInit_cgillespy(void)
{
    return PyModule_Create(&cgillespy_module_def);
}
