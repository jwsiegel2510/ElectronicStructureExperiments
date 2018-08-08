 /* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains code which tests a tangent space quadratic preconditioner on the Stiefel manifold. 
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

#include "stiefel_quadratic_preconditioner.h"
#include "stiefel_cayley_retraction.h"
#include "Eigen/Dense"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>

namespace {
using ::Eigen::MatrixXd;
using ::Eigen::IOFormat;
using ::optimization::preconditioners::StiefelQuadraticPreconditioner;
using ::optimization::retractions::StiefelCayleyRetraction;

struct DiagonalQuadratic {

void apply(MatrixXd& X) const {
	for (int j = 0; j < X.cols(); ++j) {
		for (int i = 0; i < X.rows(); ++i) {
			X(i,j) *= (i + 1);
		}
	}
}

void invert(MatrixXd& X) const {
	for (int j = 0; j < X.cols(); ++j) {
		for (int i = 0; i < X.rows(); ++i) {
			X(i,j) /= (i + 1);
		}
	}
}

};

}

int main() {
	srand(time(NULL)); // initialize random seed.
	const IOFormat fmt(-1, 1, "\t", " \n ", "(", ")", "\n", "\n");
	MatrixXd P(20, 5);
	MatrixXd V(20, 5);
	StiefelCayleyRetraction<MatrixXd, MatrixXd> retraction;
	retraction.generate_random_point(P);
	for (int j = 0; j < 5; ++j) {
		for (int i =  0; i < 20; ++i) {
                        double theta = 2 * M_PI * ((double) rand()) / (RAND_MAX);
                        double r = ((double) rand()) / (RAND_MAX);
			V(i,j) = sqrt(-log(r)) * sin(theta); 
		}
	}
	MatrixXd PreV = V;
	DiagonalQuadratic Op;
	StiefelQuadraticPreconditioner<MatrixXd, MatrixXd, DiagonalQuadratic> preconditioner(Op);
	preconditioner.precondition(PreV, P);
	if ((PreV - V).norm() < 1e-7) {
		printf("Preconditioned Point rejected, may want to run test again.\n");
	} else {
		printf("Preconditioner dual tangent space error: %lf \n", (P.transpose() * PreV + PreV.transpose() * P).norm());
		V -= (0.5 * P * (P.transpose() * V + V.transpose() * P)).eval();
		PreV += (P * (P.transpose() * PreV)).eval();
		Op.apply(PreV);
		PreV -= (0.5 * P * (P.transpose() * PreV + PreV.transpose() * P)).eval();
		printf("Preconditioner inversion error: %lf \n", (V - PreV).norm());
	}
}
