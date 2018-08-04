 /* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains code which tests the Cayley retraction. 
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

#include "stiefel_cayley_retraction.h"
#include "../../EigenLib/Eigen/Dense"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>

using ::Eigen::MatrixXd;
using ::Eigen::HouseholderQR;
using ::Eigen::IOFormat;
using ::optimization::retractions::StiefelCayleyRetraction;

int main() {
	srand(time(NULL)); // initialize random seed.
	const IOFormat fmt(2, 1, "\t", " ", "", "", "", "");
	MatrixXd P1(20, 5);
	MatrixXd V(20, 5);
	MatrixXd P2(20, 5);
	for (int i = 0; i < 20; ++i) { // initialize all points randomly
		for (int j = 0; j < 5; ++j) {
			double theta = 2 * M_PI * ((double) rand()) / (RAND_MAX);
			double r = ((double) rand()) / (RAND_MAX);
			P1(i,j) = sqrt(-log(r)) * sin(theta);
                        theta = 2 * M_PI * ((double) rand()) / (RAND_MAX);
                        r = ((double) rand()) / (RAND_MAX);
                        P2(i,j) = sqrt(-log(r)) * sin(theta);
                        theta = 2 * M_PI * ((double) rand()) / (RAND_MAX);
                        r = ((double) rand()) / (RAND_MAX);
                        V(i,j) = sqrt(-log(r)) * sin(theta);			
		}
	}
	HouseholderQR<MatrixXd> h1(P1);
	HouseholderQR<MatrixXd> h2(P2);
	P1 = h1.householderQ() * MatrixXd::Identity(P1.rows(), P1.cols());
	P2 = h2.householderQ() * MatrixXd::Identity(P2.rows(), P2.cols());
	StiefelCayleyRetraction<MatrixXd, MatrixXd> retraction;
	MatrixXd P = P1; // Test to make sure that the retraction has the correct gradient empirically.
	retraction.retract(P, V, .0001);
	MatrixXd empV = (P - P1) / .0001;
	MatrixXd trueV = V - P1 * (V.transpose() * P1);
	std::cout << (empV - trueV).format(fmt) << "\n";
	printf("Empirical Gradient Error: %lf \n", ((empV - trueV).transpose() * (empV - trueV)).trace());
	MatrixXd Q = P2;
	retraction.extrapolate(P, P1, Q, 0.0);
	printf("Extrapolation Error: %lf \n", ((P - P2).transpose() * (P - P2)).trace());
}
