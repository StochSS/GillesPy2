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

#endif

// Import solver-specific template options.
#include "template_opts.h"
