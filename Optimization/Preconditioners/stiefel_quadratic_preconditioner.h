/* This file is part of an open source library for electronic structure calculations.
 *
 *  This file contains a templated implementation of a preconditioner on the Stiefel
 *  manifold. The preconditioner inverts the Hessian of a 'separable' quadratic function defined on
 *  the tangent space of the manifold.
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

#ifndef _OPTIMIZATION_PRECONDITIONERS_STIEFEL_QUADRATIC_PRECONDITIONER__
#define _OPTIMIZATION_PRECONDITIONERS_STIEFEL_QUADRATIC_PRECONDITIONER__

#include "Eigen/Dense"
#include <cmath>

namespace optimization {
namespace preconditioners {

template<class P, class V, class Op, class EigenSolver = ::Eigen::SelfAdjointEigenSolver<MatrixXd> > 
class StiefelQuadraticPreconditioner {
private:
	const Op& operator_;
	P temp_grad;
	P temp_iterate;
	P sym_one;
	P sym_two;
	EigenSolver eigs;

public:
	StiefelQuadraticPreconditioner() = delete;
	StiefelQuadraticPreconditioner(Op& operator_in) : operator_(operator_in) {}

	double precondition(V& grad, const P& iterate) {
		temp_grad = grad;
		temp_iterate = iterate;
		operator_.invert(temp_iterate);
		sym_two = iterate.transpose() * temp_iterate();
		operator_.invert(grad);
		sym_one = iterate.transpose() * grad;
		sym_one = (sym_one + sym_one.transpose()) / 2.0;
		eigs.compute(sym_two);
		sym_one = eigs.eigenvectors().transpose() * sym_one * eigs.eigenvectors();
		auto eigenvalues = sigs.eigenvalues();
		for (int j = 0; j < sym_two.cols(); ++j) {
			for (int i = 0; i < sym_two.rows(); ++i) {
				sym_two(i,j) = -2.0 * sym_one(i,j) / (eigenvalues(i) + eigenvalues(j));
			}
		}
		grad += temp_iterate * eigs.eigenvectors() * sym_two * eigs.eigenvectors().transpose();
		
		// Only use the preconditioned direction if it is still a descent direction.
		double sq_norm = (grad.transpose() * temp_grad()).trace();
		if (sq_norm > 0) {
			grad -= 0.5 * iterate * (iterate.transpose() * grad);
			return sq_norm;
		}
		
		grad = temp_grad;
		temp_grad -= sqrt(0.5) * iterate * (grad.transpose() * iterate);
		temp_grad -= (1.0 - sqrt(0.5)) * iterate * (iterate.transpose() * grad);
		return (temp_grad.transpose() * temp_grad).trace();
	}

};

} // namespace preconditioners
} // namespace optimization

#endif
