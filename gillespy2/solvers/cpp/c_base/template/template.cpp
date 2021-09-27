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

#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include "template.h"
#include "model.h"
#include "template_definitions.h"
#include "template_defaults.h"
#include "template_params.h"

namespace Gillespy
{
	double populations[GPY_NUM_SPECIES] = GPY_INIT_POPULATIONS;
	std::vector<double> species_populations(
		populations,
		populations + sizeof(populations) / sizeof(double));

	std::string s_names[GPY_NUM_SPECIES] = 
	{
		#define SPECIES_NAME(name) #name,
		GPY_SPECIES_NAMES
		#undef SPECIES_NAME
	};

	std::vector<std::string> species_names(
		s_names,
		s_names + sizeof(s_names) / sizeof(std::string));

	int reactions[GPY_NUM_REACTIONS][GPY_NUM_SPECIES] = GPY_REACTIONS;
	std::string r_names[GPY_NUM_REACTIONS] = 
	{
		#define REACTION_NAME(name) #name,
		GPY_REACTION_NAMES
		#undef REACTION_NAME
	};

	std::vector<std::string> reaction_names(
		r_names,
		r_names + sizeof(r_names) / sizeof(std::string));

	#define VARIABLE(name, value) double name = value;
	#define CONSTANT(name, value) const double name = value;
	GPY_PARAMETER_VALUES
	#undef CONSTANT
	#undef VARIABLE

	double map_propensity(int reaction_id, const int *S)
	{
		switch (reaction_id)
		{
			#define PROPENSITY(id, func) case(id): return(func);
			GPY_PROPENSITIES
			#undef PROPENSITY

			default:
				return -1.0;
		}
	}

	double map_propensity(int reaction_id, unsigned int *S)
	{
		switch (reaction_id)
		{
			#define PROPENSITY(id, func) case(id): return(func);
			GPY_PROPENSITIES
			#undef PROPENSITY

			default:
				return -1.0;
		}
	}

	double map_propensity(int reaction_id, int *S)
	{
		switch (reaction_id)
		{
			#define PROPENSITY(id, func) case(id): return(func);
			GPY_PROPENSITIES
			#undef PROPENSITY

			default:
				return -1.0;
		}
	}

	double map_ode_propensity(int reaction_id, const std::vector<double> &S)
	{
		switch (reaction_id)
		{
			#define PROPENSITY(id, func) case(id): return(func);
			GPY_ODE_PROPENSITIES
			#undef PROPENSITY

			default:
				return -1.0;
		}
	}

	double map_ode_propensity(int reaction_id, double *S)
	{
		switch (reaction_id)
		{
			#define PROPENSITY(id, func) case(id): return(func);
			GPY_ODE_PROPENSITIES
			#undef PROPENSITY

			default:
				return -1.0;
		}
	}

	void map_variable_parameters(std::stringstream &stream)
	{
		#define VARIABLE(name, value) stream >> (name);
		#define CONSTANT(name, value)
		GPY_PARAMETER_VALUES
		#undef CONSTANT
		#undef VARIABLE
	}

	void map_variable_populations(std::stringstream &stream)
	{
		for (int spec_id = 0; spec_id < GPY_NUM_SPECIES; ++spec_id)
		{
			stream >> species_populations[spec_id];
		}
	}

	template <typename T>
	void add_reactions(Model<T> &model)
	{
		unsigned int rxn_i;
		unsigned int spec_i;

		// This is not ideal; creates an unintended side-effect!
		// Replace this with some method of "initializing" a simulation object.
		model.number_reactions = GPY_NUM_REACTIONS;
		model.number_species = GPY_NUM_SPECIES;

		for (rxn_i = 0; rxn_i < GPY_NUM_REACTIONS; ++rxn_i)
		{
			for (spec_i = 0; spec_i < GPY_NUM_SPECIES; ++spec_i)
			{
				model.reactions[rxn_i].species_change[spec_i] = reactions[rxn_i][spec_i];
			}
		}

		model.update_affected_reactions();
	}

	template void add_reactions<double>(Model<double> &model);
	template void add_reactions<unsigned int>(Model<unsigned int> &model);
	template void add_reactions<int>(Model<int> &model);
}
