#ifndef TAUHYBRIDCSOLVER_H
#define TAUHYBRIDCSOLVER_H
#include "model.h"

namespace Gillespy::TauHybrid {
    void TauHybridCSolver(Gillespy::Simulation* simulation, const double tau_tol);
}

#endif 