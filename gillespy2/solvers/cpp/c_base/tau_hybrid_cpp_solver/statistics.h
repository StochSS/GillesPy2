#pragma once

#include "HybridModel.h"
#include <vector>
#include <map>

namespace Gillespy::TauHybrid::Statistics
{

    void init_species_mode(const Model &model, Simulation &simulation);

    void partition_species(
        const Model &model,
        const std::vector<double> &propensity_values,
        std::vector<hybrid_state> current_state,
        double tau_step,
        double current_time,
        std::map<int, bool> &det_species
    );

    std::pair<std::map<std::string, int>, double> get_reactions(
        const Model *model,
        const std::vector<double> &propensity_values,
        double tau_step,
        double current_time,
        double save_time
    );

}
