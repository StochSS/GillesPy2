//
// Created by josh on 4/6/2023.
//

#pragma once

#include <functional>

namespace Gillespy
{
    template <typename PType>
    struct ModelContext
    {
    public:
        ModelContext(const std::function<double(unsigned int, PType*, double*, const double*)> &map_propensity,
                     const std::function<double(unsigned int, PType*, double*, const double*)> &map_ode_propensity);
        std::function<double(unsigned int, PType*, double*, const double*)> m_map_propensity;
        std::function<double(unsigned int, PType*, double*, const double*)> m_map_ode_propensity;
        std::function<double*(int*)> m_get_variables;
        std::function<double*(int*)> m_get_constants;
    };

    extern template struct ModelContext<unsigned int>;
}
