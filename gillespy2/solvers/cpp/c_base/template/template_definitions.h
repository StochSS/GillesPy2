#pragma once

// To set propensities, define GPY_PROPENSITY_RULES:
// g++ -D'GPY_PROPENSITY_RULES(x)=GPY_MAKE_PROPENSITY(0,x);GPY_MAKE_DEFAULT(-1)'
// Where MAKE_PROPENSITY(id,x); can be repeated as many times as necessary.

#define GPY_MAKE_PROPENSITY(id, rule) case id: return (rule)
#define GPY_MAKE_RATERULE(id, rule) case id: return (rule)
#define GPY_MAKE_DEFAULT(rule) default: return (rule)

#ifndef GPY_PROPENSITY_RULES
#define GPY_PROPENSITY_RULES(x) GPY_MAKE_DEFAULT(-2)
#endif

#ifndef GPY_RATE_RULES
#define GPY_RATE_RULES(x) GPY_MAKE_DEFAULT(-1.0)
#endif

#ifndef GPY_INIT_POPULATIONS
#define GPY_NUM_SPECIES 1
#define GPY_INIT_POPULATIONS {0}
#endif

#ifndef GPY_SPECIES_NAMES
#define GPY_SPECIES_NAMES {"x"}
#endif

#ifndef GPY_REACTION_NAMES
#define GPY_REACTION_NAMES {"r"}
#endif
