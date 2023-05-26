//
// Created by josh on 4/17/2023.
//
#pragma once

#include "model_context.h"
#include "model.h"

namespace Gillespy {

    template <typename PType>
    ModelContext<PType> *make_extension_context(const Model<PType> &model);

    extern template ModelContext<int> *make_extension_context(const Model<int> &model);
    extern template ModelContext<unsigned int> *make_extension_context(const Model<unsigned int> &model);
    extern template ModelContext<double> *make_extension_context(const Model<double> &model);
}
