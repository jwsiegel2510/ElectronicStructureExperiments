/* This file is part of an open source library for electronic structure calculations.
 *
 *  This file contains a templated implementation of a preconditioner on the Stiefel
 *  manifold. The preconditioner inverts the Hessian of a 'separable' diagonal quadratic function defined on
 *  the manifold.
 *
 *  The template parameters are classes representing a point and a tangent/dual tangent vector.
 *  Also, a class which inverts the quadratic on each column of a matrix must be provided, in
 *  addition to a class which performs a symmetric eigendecomposition. By default this is taken to
 *  be Eigen's solver since it is assumed that P and V will be eigen matrix classes. This default
 *  must be overriden if P and V are a different matrix class.
 *
 *  Points and tangent vectors are assumed to be represented by a matrix class with the same syntax
 *  as Eigen matrices.
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


#ifndef _OPTIMIZATION_PRECONDITIONERS_STIEFEL_DIAGONAL_PRECONDITIONER__
#define _OPTIMIZATION_PRECONDITIONERS_STIEFEL_DIAGONAL_PRECONDITIONER__

#endif
