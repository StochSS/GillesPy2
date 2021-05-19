#pragma once

#include <vector>
#include <string>

#include "../model.h"
#include "XTemplateDefinitions.h"

unsigned int populations[GPY_NUM_SPECIES] = GPY_INIT_POPULATIONS;
std::string s_names[GPY_NUM_SPECIES] = GPY_SPECIES_NAMES;

int reactions[GPY_NUM_REACTIONS][GPY_NUM_SPECIES] = GPY_REACTIONS;
std::string r_names[GPY_NUM_REACTIONS] = GPY_REACTION_NAMES;

namespace Gillespy {
    inline double map_propensity(int reaction_id, std::vector<unsigned int> state) {
        switch (reaction_id) {
            #define PROPENSITY(id, func) case(id): return(func)
            GPY_PROPENSITIES
            #undef PROPENSITY(id, func)
        }
    }

    inline double map_ode_propensity(int reaction_number, const std::vector<double> state)
    {
        switch (reaction_number) {
            #define RATE(id, func)
            GPY_RATES
            #undef RATE(id, func)
        }
    }

    inline void add_reactions(Model &model)
    {
        unsigned int rxn_i;
        unsigned int spec_i;

        for (rxn_i = 0; rxn_i < GPY_NUM_REACTIONS; ++rxn_i) {
            for (spec_i = 0; spec_i < GPY_NUM_SPECIES; ++spec_i) {
                model.reactions[rxn_i].species_change[spec_i] = reactions[rxn_i][spec_i];
            }
        }
        
        model.update_affected_reactions();
    }
}