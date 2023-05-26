/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2023 GillesPy2 developers.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include "template.h"
#include "model.h"

namespace Gillespy
{
    std::vector<double> species_populations = { 100.0f };
    std::vector<std::string> species_names = { "x" };
    std::vector<std::string> reaction_names = { "Rx" };

    template <typename PType>
    double map_ssa_propensity(unsigned int reaction_id, PType *S, double *P, const double *C)
    {
        return 0.2 * S[0];
    }

    template <typename PType>
    double map_ode_propensity(unsigned int reaction_id, PType *S, double *P, const double *C)
    {
        return 0.2 * S[0];
    }

    double *get_variables(int *num_variables)
    {
        *num_variables = 0;
        return new double[0];
    }

    double *get_constants(int *num_constants)
    {
        *num_constants = 0;
        return new double[0];
    }

    template <typename T>
    void add_reactions(Model<T> &model)
    {
        model.number_reactions = 1;
        model.number_species = 1;
        model.reactions[0].species_change[0] = -1;
        model.reactions[0].reactants_change[0] = 1;
        model.reactions[0].products_change[0] = 0;
    }

    void map_variable_parameters(std::stringstream &stream)
    {
        // do nothing
    }

    void map_variable_populations(std::stringstream &stream)
    {
        // do nothing
    }

    template double map_ssa_propensity(unsigned int reaction_id, unsigned int *S, double *P, const double *C);
    template double map_ode_propensity(unsigned int reaction_id, unsigned int *S, double *P, const double *C);
    template void add_reactions<unsigned int>(Model<unsigned int> &model);
}
