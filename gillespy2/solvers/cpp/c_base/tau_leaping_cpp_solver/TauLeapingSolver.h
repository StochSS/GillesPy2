#pragma once
#include "model.h"

namespace Gillespy {
	void tau_leaper(Gillespy::Simulation<unsigned int> *simulation, const double tau_tol);
}
