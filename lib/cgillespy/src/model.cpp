/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2023 GillesPy2 developers.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "model.h"
#include <algorithm>

namespace Gillespy
{
    template <typename PType>
    Model<PType>::Model(
            const ModelContext<PType> &context,
            std::vector<std::string> species_names,
            std::vector<double> species_populations,
            std::vector<std::string> reaction_names)
        : number_species(species_names.size()),
          number_reactions(reaction_names.size())
    {
        species = std::make_unique<Species<PType>[]>(number_species);
        reactions = std::make_unique<Reaction<PType>[]>(number_reactions);

        int num_constants;
        double *constants = context.m_get_constants(&num_constants);
        for (int i = 0; i < num_constants; ++i) {
            m_initial_constants.push_back(constants[i]);
        }

        int num_variables;
        double *variables = context.m_get_variables(&num_variables);
        for (int i = 0; i < num_variables; ++i)
        {
            m_initial_variables.push_back(variables[i]);
        }

        for (unsigned int i = 0; i < number_species; i++) {
            species[i].id = i;
            species[i].initial_population = static_cast<PType>(species_populations[i]);
            species[i].name = species_names[i];
        }

        for (unsigned int reaction = 0; reaction < number_reactions; reaction++) {
            reactions[reaction].name = reaction_names[reaction];
            reactions[reaction].species_change = std::make_unique<int[]>(number_species);
            reactions[reaction].reactants_change = std::make_unique<int[]>(number_species);
            reactions[reaction].products_change = std::make_unique<int[]>(number_species);
            reactions[reaction].m_map_ssa_propensity = context.m_map_propensity;
            reactions[reaction].m_map_ode_propensity = context.m_map_ode_propensity;

            for (unsigned int species = 0; species < number_species; species++) {
                reactions[reaction].species_change[species] = 0;
                reactions[reaction].reactants_change[species] = 0;
                reactions[reaction].products_change[species] = 0;
            }

            reactions[reaction].affected_reactions = std::vector<unsigned int>();
        }
    }

    template <typename PType>
    std::vector<double> Model<PType>::copy_variables() const
    {
        return m_initial_variables;
    }

    template <typename PType>
    std::vector<double> Model<PType>::copy_constants() const
    {
        return m_initial_constants;
    }

    template <typename PType>
    void Model<PType>::update_affected_reactions() {
        // Clear affected_reactions for each reaction.
        for (unsigned int i = 0; i < number_reactions; i++) {
            reactions[i].affected_reactions.clear();
        }

        // Check all reactions for commong species changes -> affected reactions.
        for (unsigned int r1 = 0; r1 < number_reactions; r1++) {
            for (unsigned int r2 = 0; r2 < number_reactions; r2++) {
                for (unsigned int s = 0; s < number_species; s++) {
                     if(reactions[r1].species_change[s]  != 0  &&
                        reactions[r2].reactants_change[s] > 0 ){
                        reactions[r1].affected_reactions.push_back(r2);
                    }
                }
            }
        }
    }

    LogStream &Logger::info()
    {
        switch (m_log_level)
        {
        case LogLevel::INFO:
            return m_stderr;
        default:
            return m_null;
        }
    }

    LogStream &Logger::warn()
    {
        switch (m_log_level)
        {
        case LogLevel::WARN:
        case LogLevel::INFO:
            return m_stderr;
        default:
            return m_null;
        }
    }

    LogStream &Logger::err()
    {
        switch (m_log_level)
        {
        case LogLevel::CRIT:
        case LogLevel::SILENT:
            return m_null;
        default:
            return m_stderr;
        }
    }

    LogStream &Logger::crit()
    {
        switch (m_log_level)
        {
        case LogLevel::SILENT:
            return m_null;
        default:
            return m_stderr;
        }
    }

    // DETERMINISTIC SIMULATIONS: explicit instantation of real-valued data structures.
    template struct Species<double>;
    template struct Reaction<double>;
    template class Model<double>;

    // STOCHASTIC SIMULATIONS: explicit instantiation of discrete-valued data structures.
    template struct Species<unsigned int>;
    template struct Species<int>;
    template struct Reaction<unsigned int>;
    template struct Reaction<int>;
    template class Model<unsigned int>;
    template class Model<int>;
}
