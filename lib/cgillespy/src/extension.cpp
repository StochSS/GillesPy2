//
// Created by josh on 2/28/2023.
//

#include "model.h"
#include "simulation.h"
#include "extension.h"
#include "ssa_cpp_solver/SSASolver.h"

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
            (char*) "species", // { [species_name]: <initial_value> }
            (char*) "reactions", // { [reaction_name]: { reactants: { [species]: <value> }, products: { [species]: <value> } }
            (char*) "constants", // { [constant_name]: <value> }
            (char*) "variables", // { [variable_name]: <value> }
            NULL
    };

    PyObject *py_species_dict;
    PyObject *py_reactions_dict;
    PyObject *py_constants_list;
    PyObject *py_variables_list;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOO", keywords, &py_species_dict, &py_reactions_dict, &py_constants_list, &py_variables_list))
    {
        return PyUnicode_FromString("Failed to parse arguments");
    }

    // Validate that each argument conforms to the expected API
    if (!(PyDict_Check(py_species_dict) &&
          PyDict_Check(py_reactions_dict) &&
          PyDict_Check(py_constants_list) &&
          PyDict_Check(py_variables_list)
       ))
    {
        return PyUnicode_FromString("Incorrect usage");
    }

    PyObject *key;
    PyObject *value;
    Py_ssize_t position = 0;
    std::vector<std::string> species_names;
    std::vector<double> initial_populations;
    std::vector<std::string> reaction_names;

    while (PyDict_Next(py_species_dict, &position, &key, &value))
    {
        // Verify that the key is a string, and the value is a number
        if (!PyUnicode_Check(key) || !PyFloat_Check(value))
        {
            return PyUnicode_FromString("Invalid key/value in species dict");
        }

        const char *species_name = PyUnicode_AsUTF8(key);
        double population_value = PyFloat_AsDouble(value);
        std::cout << std::string(species_name) << " = " << population_value << std::endl;
        species_names.emplace_back(species_name);
        initial_populations.emplace_back(population_value);
    }

    position = 0;
    while (PyDict_Next(py_reactions_dict, &position, &key, &value))
    {
        if (!PyUnicode_Check(key) || !PyDict_Check(value))
        {
            return PyUnicode_FromString("Invalid key/value in reactions dict");
        }

        const char *reaction_name = PyUnicode_AsUTF8(key);
        std::cout << std::string(reaction_name) << " (rxn)" << std::endl;
        reaction_names.emplace_back(reaction_name);
    }


    position = 0;
    std::vector<std::vector<int>> react_stoch;
    std::vector<std::vector<int>> prod_stoch;
    std::vector<std::function<double(unsigned int*, double*, const double*)>> propensity_functions;
    while (PyDict_Next(py_reactions_dict, &position, &key, &value))
    {
        std::vector<size_t> dependent_species;
        std::vector<int> react;
        std::vector<int> prod;
        PyObject *py_reactants_dict = PyDict_GetItemString(value, "reactants");
        PyObject *py_products_dict = PyDict_GetItemString(value, "products");
        PyObject *py_rate_index = PyDict_GetItemString(value, "rate_index");
        if (py_reactants_dict == nullptr || py_products_dict == nullptr || py_rate_index == nullptr)
        {
            return PyUnicode_FromString("reactants/products/index are null");
        }
        if (!PyList_Check(py_reactants_dict) || !PyList_Check(py_products_dict) || !PyLong_Check(py_rate_index))
        {
            Py_DECREF(py_reactions_dict);
            Py_DECREF(py_products_dict);
            return PyUnicode_FromString("reactants/products are not iterable or index is not int");
        }

        Py_ssize_t index;
        for (index = 0; index < PyList_Size(py_reactants_dict); ++index)
        {
            PyObject *py_reaction_item = PyList_GetItem(py_reactants_dict, index);
            if (!PyLong_Check(py_reaction_item))
            {
                Py_DECREF(py_reaction_item);
                return PyUnicode_FromString("reactant value not an int");
            }
            long stoch_value = PyLong_AsLong(py_reaction_item);
            if (stoch_value != 0) {
                dependent_species.push_back(index);
            }
            react.emplace_back(stoch_value);

            Py_DECREF(py_reaction_item);
        }

        for (index = 0; index < PyList_Size(py_products_dict); ++index)
        {
            PyObject *py_reaction_item = PyList_GetItem(py_products_dict, index);
            if (!PyLong_Check(py_reaction_item))
            {
                Py_DECREF(py_reaction_item);
                return PyUnicode_FromString("product value not an int");
            }
            long stoch_value = PyLong_AsLong(py_reaction_item);
            prod.emplace_back(stoch_value);

            Py_DECREF(py_reaction_item);
        }

        long rate_index = PyLong_AsLong(py_rate_index);
        propensity_functions.emplace_back([dependent_species, rate_index](unsigned int *current_state, double *variables, const double *constants) -> double {
            double result = variables[rate_index];
            for (const auto spec_i : dependent_species)
            {
                result *= current_state[spec_i];
            }
            return result;
        });

        react_stoch.emplace_back(react);
        prod_stoch.emplace_back(prod);
        Py_DECREF(py_reactions_dict);
        Py_DECREF(py_products_dict);
    }

    std::vector<double> base_variables;
    std::vector<double> base_constants;
    position = 0;
    while (PyDict_Next(py_variables_list, &position, &key, &value))
    {
        if (!PyFloat_Check(value))
        {
            return PyUnicode_FromString("bad variable value");
        }

        double variable_value = PyFloat_AsDouble(value);
        base_variables.emplace_back(variable_value);
    }

    position = 0;
    while (PyDict_Next(py_constants_list, &position, &key, &value))
    {
        if (!PyFloat_Check(value))
        {
            return PyUnicode_FromString("bad constant value");
        }

        double constant_value = PyFloat_AsDouble(value);
        base_constants.emplace_back(constant_value);
    }

    std::function<double(unsigned int, unsigned int*, double*, const double*)> map_propensity = [&propensity_functions](unsigned int reaction_id, unsigned int *current_state, double *variables, const double *constants) -> double {
        double result = propensity_functions[reaction_id](current_state, variables, constants);
        return result;
    };
    Gillespy::ModelContext<unsigned int> model_context(map_propensity, map_propensity);
    model_context.m_get_variables = [&base_variables](int *num_variables) -> double* {
        *num_variables = base_variables.size();
        return base_variables.data();
    };
    model_context.m_get_constants = [&base_constants](int *num_constants) -> double* {
        *num_constants = base_constants.size();
        return base_constants.data();
    };
    Gillespy::Model<unsigned int> model(model_context, species_names, initial_populations, reaction_names);
    Gillespy::Simulation<unsigned int> simulation;
//    simulation.output_interval = 20;
    simulation.current_time = 0;
    simulation.model = &model;
    simulation.end_time = 20;
    simulation.random_seed = time(nullptr);
    simulation.number_timesteps = 21;
    simulation.number_trajectories = 1;
    std::cout << "DBG rxn# = " << model.number_reactions << std::endl;
    Gillespy::init_simulation(&model, simulation);
    for (int rxn_i = 0; rxn_i < model.number_reactions; ++rxn_i)
    {
        model.reactions[rxn_i].id = rxn_i;
        for (int spec_i = 0; spec_i < model.number_species; ++spec_i)
        {
            model.reactions[rxn_i].reactants_change[spec_i] = react_stoch[rxn_i][spec_i];
            model.reactions[rxn_i].products_change[spec_i] = prod_stoch[rxn_i][spec_i];
            model.reactions[rxn_i].species_change[spec_i] = model.reactions[rxn_i].products_change[spec_i] - model.reactions[rxn_i].reactants_change[spec_i];
        }
    }
    model.update_affected_reactions();
    Gillespy::ssa_direct(&simulation);
    simulation.output_buffer_final(std::cout);

    return PyLong_FromLong(0);
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
