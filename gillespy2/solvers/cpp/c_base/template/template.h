/* TEMPLATE.H
 * 
 * This is the header file which defines the interface of the template.
 * Includes functions for loading and defining simulation parameters.
 */

#include "model.h"

namespace Gillespy
{
    extern std::vector<unsigned int> species_populations;
    extern std::vector<std::string> species_names;

    extern int *reactions[];
    extern std::vector<std::string> reaction_names;

    extern inline double map_propensity(int reaction_id, std::vector<unsigned int> state);
    extern inline double map_ode_propensity(int reaction_id, std::vector<double> state);
    extern inline void add_reactions(Model &model);

}
