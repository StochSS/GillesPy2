//
// Created by josh on 4/6/2023.
//

#include "model_context.h"

template <typename PType>
Gillespy::ModelContext<PType>::ModelContext(const std::function<double(unsigned int, PType*, double*, const double*)> &map_propensity,
                                            const std::function<double(unsigned int, PType*, double*, const double*)> &map_ode_propensity)
    : m_map_propensity(map_propensity),
      m_map_ode_propensity(map_ode_propensity)
{

}

template struct Gillespy::ModelContext<unsigned int>;
