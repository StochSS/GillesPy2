#pragma once

#include "template_definitions.h"
#include "../model.h"
#include <vector>
#include <string>

unsigned int populations[GPY_NUM_SPECIES] = GPY_INIT_POPULATIONS;
std::string s_names[GPY_NUM_SPECIES] = GPY_SPECIES_NAMES;

int reactions[GPY_NUM_REACTIONS][GPY_NUM_SPECIES] = GPY_REACTIONS;
std::string r_names[GPY_NUM_REACTIONS] = GPY_REACTION_NAMES;

namespace Gillespy {

    inline double map_propensity(int reaction_number, std::vector<unsigned int> state)
    {
        switch (reaction_number) {
            GPY_PROPENSITY_RULES(reaction_number);
        }
    }

    inline double map_ode_propensity(int reaction_number, const std::vector<double> state)
    {
        switch (reaction_number) {
            GPY_RATE_RULES(reaction_number);
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