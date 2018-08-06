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
#include <cstdio>
#include <iostream>

namespace {
	using ::Eigen::MatrixXd;
	using ::Eigen::IOFormat;
	using ::optimization::methods::accelerated_gradient_descent;
	using ::optimization::retractions::StiefelCayleyRetraction;
	using ::compressed_modes::CompressedModesObjective;

	struct LaplaceOperator { // 1D Laplacian operator on Eigen matrices.
		MatrixXd apply(const MatrixXd& input) const {
			MatrixXd output(input);
			if (input.rows() <= 1) return 2.0 * output;
			for (int j = 0; j < input.cols(); ++j) {
				for (int i = 0; i < input.rows(); ++i) {
					if (i == 0) {
						output(i,j) = 2 * input(i,j) - input(i+1,j);
					} else if (i == input.rows() - 1) {
						output(i,j) = 2 * input(i,j) - input(i-1,j);
					} else {
						output(i,j) = 2 * input(i,j) - input(i+1,j) - input(i-1,j);
					}
				}
			}
			return output;
		}
	};	
}

int main() {
	double mu = .1;
	double epsilon = .001;
	const IOFormat fmt(-1, 1, "\t", " \n ", "(", ")", "\n", "\n");	
	MatrixXd iterate(50,5);
	LaplaceOperator op_;
	CompressedModesObjective<MatrixXd, MatrixXd, LaplaceOperator> objective(op_, mu, epsilon);
	StiefelCayleyRetraction<MatrixXd, MatrixXd> retraction;
	retraction.generate_random_point(iterate);
	accelerated_gradient_descent(iterate, objective, retraction);
	std::cout << iterate.format(fmt) << "\n";
}
