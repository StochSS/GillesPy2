/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2022 GillesPy2 developers.
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

/* TEMPLATE.H
 *
 * This is the header file which defines the interface of the template.
 * Includes functions for loading and defining simulation parameters.
 */

#include <sstream>
#include <vector>

namespace Gillespy
{
    template <typename PType>
    struct Model;

    extern std::vector<double> species_populations;
    extern std::vector<std::string> species_names;
    extern std::vector<std::string> reaction_names;

    double map_propensity(unsigned int reaction_id, int *state, double *parameters, const double *constants);
    double map_propensity(unsigned int reaction_id, unsigned int *S, double *parameters, const double *constants);
    double map_propensity(unsigned int reaction_id, double *S, double *parameters, const double *constants);

    double map_ssa_propensity(unsigned int reaction_id, int *state, double *parameters, const double *constants);
    double map_ssa_propensity(unsigned int reaction_id, unsigned int *state, double *parameters, const double *constants);
    double map_ssa_propensity(unsigned int reaction_id, double *state, double *parameters, const double *constants);
    double map_ode_propensity(unsigned int reaction_id, int *state, double *parameters, const double *constants);
    double map_ode_propensity(unsigned int reaction_id, unsigned int *state, double *parameters, const double *constants);
    double map_ode_propensity(unsigned int reaction_id, double *state, double *parameters, const double *constants);

    template <typename T>
    void add_reactions(Model<T> &model);

    double *get_variables(int *num_variables);
    double *get_constants(int *num_constants);

    void map_variable_parameters(std::stringstream &stream);
    void map_variable_populations(std::stringstream &stream);

    extern template void add_reactions<double>(Model<double> &model);
    extern template void add_reactions<unsigned int>(Model<unsigned int> &model);
    extern template void add_reactions<int>(Model<int> &model);
}
