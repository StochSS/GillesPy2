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

#include "hybrid_template.h"
#include "template_params.h"

namespace Gillespy
{
	namespace TauHybrid
	{
		void map_species_modes(std::vector<HybridSpecies> &species)
		{
			#define SPECIES_MODE(spec_id, user_min, spec_mode, boundary_mode) \
			species[spec_id].user_mode = spec_mode; \
			species[spec_id].switch_min = user_min; \
			species[spec_id].boundary_condition = boundary_mode;
			#define CONTINUOUS_MODE SimulationState::CONTINUOUS
			#define DISCRETE_MODE   SimulationState::DISCRETE
			#define DYNAMIC_MODE    SimulationState::DYNAMIC
			#define BOUNDARY true
			#define STANDARD false
			GPY_HYBRID_SPECIES_MODES
			#undef STANDARD
			#undef BOUNDARY
			#undef DYNAMIC_MODE
			#undef DISCRETE_MODE
			#undef CONTINUOUS_MODE
			#undef SPECIES_MODE
		}

		void map_rate_rules(std::vector<HybridSpecies> &species)
		{
			#define RATE_RULE(spec_id, rate_rule) species[spec_id].diff_equation.rate_rules.push_back([](double t, double *S) { return (rate_rule); });
			GPY_RATE_RULES
			#undef RATE_RULE
		}

	}
}
