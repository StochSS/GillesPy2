#include "hybrid_template.h"

namespace Gillespy::TauHybrid
{
	void map_species_modes(std::vector<HybridSpecies> &species)
	{
		#define SPECIES_MODE(spec_id, spec_mode) species[spec_id].user_mode = spec_mode;
		#define CONTINUOUS SimulationState::CONTINUOUS
		#define DISCRETE   SimulationState::DISCRETE
		#define DYNAMIC    SimulationState::DYNAMIC
		GPY_HYBRID_SPECIES_MODES
		#undef DYNAMIC
		#undef DISCRETE
		#undef CONTINUOUS
		#undef SPECIES_MODE
	}
}
