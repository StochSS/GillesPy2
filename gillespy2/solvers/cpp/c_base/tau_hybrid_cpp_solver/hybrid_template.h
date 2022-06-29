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

#pragma once

#define GPY_SOLVER_HYBRID

#include "template.h"
#include "HybridModel.h"
#include "template_defaults.h"

#include <functional>
#include <initializer_list>

#ifndef GPY_HYBRID_EVENTS
#define GPY_HYBRID_EVENTS
#define GPY_HYBRID_NUM_EVENTS 0
#endif

#ifndef GPY_HYBRID_EVENT_ASSIGNMENTS
#define GPY_HYBRID_EVENT_ASSIGNMENTS
#define GPY_HYBRID_NUM_EVENT_ASSIGNMENTS 0
#endif

namespace Gillespy
{
    namespace TauHybrid
    {

        void map_species_modes(std::vector<HybridSpecies> &species);
        void map_rate_rules(std::vector<HybridSpecies> &species);

    }
}
