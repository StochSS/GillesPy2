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

namespace Gillespy::TauHybrid
{
	void map_species_modes(std::vector<HybridSpecies> &species)
	{
		#define SPECIES_MODE(spec_id, spec_mode, user_min) \
		species[spec_id].user_mode = spec_mode; \
		species[spec_id].switch_min = user_min;
		#define CONTINUOUS_MODE SimulationState::CONTINUOUS
		#define DISCRETE_MODE   SimulationState::DISCRETE
		#define DYNAMIC_MODE    SimulationState::DYNAMIC
		GPY_HYBRID_SPECIES_MODES
		#undef DYNAMIC_MODE
		#undef DISCRETE_MODE
		#undef CONTINUOUS_MODE
		#undef SPECIES_MODE
	}
}
