#pragma once

#ifndef GPY_PROPENSITIES
#define GPY_PROPENSITIES PROPENSITY(0, x)
#endif

#ifndef GPY_RATES
#define GPU_RATES RATE(0, x)
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
#define GPY_SPECIES_NAMES {"x"}
#endif

#ifndef GPY_REACTION_NAMES
#define GPY_REACTION_NAMES {"r"}
#endif
