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
