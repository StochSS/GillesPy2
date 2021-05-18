#pragma once

#include "template_definitions.h"
#include <vector>
#include <string>

unsigned int populations[] = GPY_INIT_POPULATIONS;
std::string s_names[] = GPY_SPECIES_NAMES;
std::string r_names[] = GPY_REACTION_NAMES;

inline double map_propensity(int reaction_number, std::vector<unsigned int> state)
{
    switch (reaction_number) {
        GPY_PROPENSITY_RULES(reaction_number);
    }
}

inline double map_rate(int reaction_number, double t, std::vector<double> state)
{
    switch (reaction_number) {
        GPY_RATE_RULES(reaction_number);
    }
}
