#ifndef TAU_LEAPER_H
#define TAU_LEAPER_H
#include "model.h"

namespace Gillespy{
    void tau_leaper(Gillespy::Simulation<unsigned int> *simulation, const double tau_tol);
}

#endif // TAU_LEAPER_H
