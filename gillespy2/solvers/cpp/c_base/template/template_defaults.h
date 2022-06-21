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

// To set propensities: 
// g++ -D GPY_PROPENSITIES='PROPENSITY(0, x) PROPENSITY(1, x + 1) PROPENSITY(2, x + 2)'

// To set rate rules:
// g++ -D GPY_RATES='RATE(0, x) RATE(1, x + 1) RATE(2, x + 2)'

#ifndef GPY_PROPENSITIES
#define GPY_PROPENSITIES
#endif

#ifndef GPY_ODE_PROPENSITIES
#define GPY_ODE_PROPENSITIES GPY_PROPENSITIES
#endif

#ifndef GPY_PARAMETER_VALUES
#define GPY_PARAMETER_VALUES
#define GPY_PARAMETER_NUM_VARIABLES 0
#define GPY_PARAMETER_NUM_CONSTANTS 0
#endif

#ifndef GPY_INIT_POPULATIONS
#define GPY_NUM_SPECIES 1
#define GPY_INIT_POPULATIONS {0}
#endif

#ifndef GPY_REACTIONS
#define GPY_NUM_REACTIONS 1
#define GPY_REACTIONS { {0} }
#endif

#ifndef GPY_VOLUME
#define GPY_VOLUME 1.0
#endif

#ifndef GPY_SPECIES_NAMES
#define GPY_SPECIES_NAMES
#endif

#ifndef GPY_REACTION_NAMES
#define GPY_REACTION_NAMES
#endif

// ===============================================================
// ================ HYBRID SOLVER OPTION DEFAULTS ================
// ===============================================================

/* GPY_HYBRID_SPECIES_MODES: Default, user-provided flags for how each species is to be represented.
 * Populate each SPECIES_MODE() with two arguments: species ID and species mode.
 * Possible values for species mode are: DISCRETE / CONTINUOUS / DYNAMIC
 * 
 * 
 * #define GPY_HYBRID_SPECIES_MODES \
 * SPECIES_MODE(0, DISCRETE_MODE) \
 * SPECIES_MODE(1, CONTINUOUS_MODE) \
 * SPECIES_MODE(2, DYNAMIC_MODE)
 */
#ifdef GPY_SOLVER_HYBRID

#ifndef GPY_HYBRID_SPECIES_MODES
#define GPY_HYBRID_SPECIES_MODES
#endif

#ifndef GPY_RATE_RULES
#define GPY_RATE_RULES
#endif

#endif

// Import solver-specific template options.
#include "template_opts.h"
