/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2021 GillesPy2 developers.
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

#include "model.h"

namespace Gillespy
{
	extern std::vector<double> species_populations;
	extern std::vector<std::string> species_names;
	extern std::vector<std::string> reaction_names;

	double map_propensity(int reaction_id, const int *state);
	double map_propensity(int reaction_id, unsigned int *S);
	double map_propensity(int reaction_id, int *S);
	double map_ode_propensity(int reaction_id, const std::vector<double> &state);
	double map_ode_propensity(int reaction_id, double *S);

	template <typename T>
	void add_reactions(Model<T> &model);

	void map_variable_parameters(std::stringstream &stream);
	void map_variable_populations(std::stringstream &stream);

	extern template void add_reactions<double>(Model<double> &model);
	extern template void add_reactions<unsigned int>(Model<unsigned int> &model);
	extern template void add_reactions<int>(Model<int> &model);
}
