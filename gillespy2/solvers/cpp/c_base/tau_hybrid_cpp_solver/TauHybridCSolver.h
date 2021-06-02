#ifndef TAUHYBRIDCSOLVER_H
#define TAUHYBRIDCSOLVER_H
#include "HybridModel.h"

namespace Gillespy::TauHybrid {
    void TauHybridCSolver(HybridSimulation* simulation, const double tau_tol);
}

#endif 