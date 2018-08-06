/* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains code which calculates compressed modes for the laplacian in 1D, based on the 
 *  parameters given in an XML file passed to the program.
 *
 *  Copyright (C) 2018 Jonathan W. Siegel
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../XML-ParameterList/XML_ParameterListArray.h"
#include "../EigenLib/Eigen/Dense"
#include "../Optimization/Methods/accelerated_gradient_descent.h"
#include "../Optimization/Retractions/stiefel_cayley_retraction.h"
#include "compressed_modes_objective.h"

namespace {
	
}
