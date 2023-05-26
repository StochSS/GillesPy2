//
// Created by josh on 4/17/2023.
//

#include "template_extension.h"
#include <vector>
#include <functional>

template <typename PType>
Gillespy::ModelContext<PType> *Gillespy::make_extension_context(const Gillespy::Model<PType> &model)
{
    std::vector<std::function<double(PType*, double*, double*)>> propensity_functions;

    for (Gillespy::Reaction<PType> &reaction : model.reactions.get())
    {
        std::vector<size_t> dependent_species;

        // go through list of reactants, mark their id
        // then add to the list of propensities
        for (int spec_i = 0; spec_i < model.number_species; ++spec_i)
        {
            int spec_dx = reaction.reactants_change[spec_i];
            if (spec_dx != 0)
                dependent_species.emplace_back(spec_i);
        }

        std::function<double(PType*, double*, double*)> propensity_function = [&dependent_species](
                PType *current_state, double *variables, double *constants) -> double
        {
            double result = 0.0;
            for (size_t dependent_spec : dependent_species)
            {
                result *= dependent_spec;
            }
        };
    }
}
