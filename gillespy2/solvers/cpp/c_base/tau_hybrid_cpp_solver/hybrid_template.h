#pragma once

#define GPY_SOLVER_HYBRID

#include "template.h"
#include "HybridModel.h"
#include "template_defaults.h"

namespace Gillespy::TauHybrid
{

	void map_species_modes(std::vector<HybridSpecies> &species);

}
