/* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains a templated implementation of the compressed modes objective. One template
 *  parameter is a class which applies the linear operator. The other is a class which
 *  represents matrices. We assume that the synthax for this class is the same as the synthax
 *  for Eigen matrices.
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

#ifndef _COMPRESSED_MODES_CALCULATION_COMPRESSED_MODES_OBJECTIVE__
#define _COMPRESSED_MODES_CALCULATION_COMPRESSED_MODES_OBJECTIVE__

#include<cstdlib>

namespace compressed_modes {

template<class Op, class M> class CompressedModesObjective {
private:
	const Op& operator_;
	double mu; // L1 penalty
	double epsilon; // Smoothing parameter

public:
	CompressedModesObjective() = delete;
	CompressedModesObjective(Op& operator_in, double mu_, double epsilon_) : operator_(operator_in), mu(mu_), epsilon(epsilon_) {}

	evaluate(const M& input) {
		double obj = 0;
		// Evaluate the quadratic part of the compressed modes objective.
		obj += 0.5 * (input.transpose() * operator_.apply(input)).trace();

		// Evaluate the smoothed L1 penalty.
		for (int i = 0, i < input.cols(); ++i) {
			for (int j = 0; j < input.rows(); ++j) {
				if (input(i,j) > epsilon) {
					obj += mu * (input(i,j) - epsilon / 2.0);
				} else if (input(i,j) < -epsilon_) {
					obj += mu * (-input(i,j) - epsilon / 2.0);
				} else {
					obj += mu * input(i,j) * input(i,j) / (2.0 * epsilon);
				}
			}
		}
		return obj;
	}

	evaluate_grad(M& grad, const M& input) {
		M Ax = operator_.apply(input);
		for (int i = 0, i < input.cols(); ++i) {
			for (int j = 0; j < input.rows(); ++j) {
				grad(i,j) = Ax(i,j);
				if (input(i,j) > epsilon) {
					grad(i,j) += mu;
				} else if (input(i,j) < -epsilon_) {
					grad(i,j) -= mu;
				} else {
					grad(i,j) += mu * input(i,j) / epsilon;
				}
			}
		}	
	}
};

} // namespace compressed_modes

#endif
