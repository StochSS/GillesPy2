#pragma once

#include "cvode.h"

namespace Gillespy::TauHybrid
{
	int rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data);
}
