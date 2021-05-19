#include <vector>
#include <string>

#include "model.h"
#include "XTemplateDefinitions.h"
#include "XTemplateDefaults.h"

namespace Gillespy {

    unsigned int populations[GPY_NUM_SPECIES] = GPY_INIT_POPULATIONS;
    std::vector<unsigned int> species_populations(
        populations,
        populations + sizeof(populations) / sizeof(unsigned int));

    std::string s_names[GPY_NUM_SPECIES] = GPY_SPECIES_NAMES;
    std::vector<std::string> species_names(
        s_names,
        s_names + sizeof(s_names) / sizeof(std::string));

    int reactions[GPY_NUM_REACTIONS][GPY_NUM_SPECIES] = GPY_REACTIONS;
    std::string r_names[GPY_NUM_REACTIONS] = GPY_REACTION_NAMES;
    std::vector<std::string> reaction_names(
        r_names,
        r_names + sizeof(r_names) / sizeof(std::string));

    double map_propensity(int reaction_id, std::vector<unsigned int> S) {
        switch (reaction_id) {
            #define PROPENSITY(id, func) case(id): return(func);
            GPY_ODE_PROPENSITIES
            #undef PROPENSITY

            default:
                return -1.0;
        }
    }

    double map_ode_propensity(int reaction_id, const std::vector<double> S)
    {
        switch (reaction_id) {
            #define PROPENSITY(id, func) case(id): return(func);
            GPY_PROPENSITIES
            #undef PROPENSITY

            default:
                return -1.0;
        }
    }

    void add_reactions(Model &model)
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