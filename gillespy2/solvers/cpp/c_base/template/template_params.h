/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2022 GillesPy2 developers.
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

/* template_params.h
 * 
 * Include this header file for access to the model-level parameter values.
 * Parameters are provided by iterating over the GPY_PARAMETER_VALUES macro.
 * The definitions for all parameters are provided in template.cpp.
 * 
 * NOTE: as with all other template files, this header is only valid
 *   when compiled along with the model, and as such CANNOT be precompiled!
 */

#pragma once

#include "template_definitions.h"
#include <cmath>

// "factorial()" function is not defined in cmath
// Can be emulated using std::tgamma
//   std::gamma(n) = factorial(n-1),
//     and so:
//   std::gamma(n+1) = factorial((n+1)-1) = factorial(n)
#ifndef factorial
#define factorial(n) std::tgamma(static_cast<int>(n)+1)
#endif

using namespace std;
