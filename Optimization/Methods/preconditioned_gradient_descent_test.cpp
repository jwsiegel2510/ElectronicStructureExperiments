/* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains code which tests the preconditioned gradient descent method on a 
 *  simple quadratic function. 
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

#include<vector>
#include<cstdlib>
#include<cstdio>
#include "test_lib.h"
#include "preconditioned_gradient_descent.h"

using std::vector;

int main() {
	int n = 10;
	vector<double> iterate(10, 1.0);
	int iteration_count = preconditioned_gradient_descent(iterate, Objective<>(), Retraction<>(), Preconditioner<>());
	if (iterate.size() != 10) {
		printf("Test Failed: Vector changed size unexpectedly.\n");
		return -1;
	}
	for (int i = 0; i < 10; ++i) {
		if (abs(iterate[i]) > 1e-4) {
			printf("Test Failed: Failure to converge.\n");
			return -1;
		}
	}
	printf("Test passed. Converged in %d iterations.\n", iteration_count);
}
