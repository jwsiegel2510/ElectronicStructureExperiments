 /* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains code which tests the gradient descent method on a simple quadratic function. 
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
#include "gradient_descent.h"

namespace {
	using std::vector;
	
	template<class P = vector<double>, class V = vector<double> > struct Objective {
		
		double evaluate(const P& input) {
			double output = 0;
			for (int i = 0; i < input.size(); ++i) {
				output += .5 * (i+1) * input[i] * input[i];
			}
			return output;
		}

		void evaluate_grad(V& grad, const P& input) {
			grad.resize(input.size());
			for (int i = 0; i < input.size(); ++i) {
				grad[i] = (i+1) * input[i];
			}
		}
	};

	template<class P = vector<double>, class V = vector<double> > struct Retraction {
		
		void retract(P& input, const V& direction, double distance) {
			for (int i = 0; i < input.size(); ++i) {
				input[i] += direction[i] * distance;
			}
		}

		double norm_sq(const V& vect, const P& point) {
			double norm_sq = 0;
			for (int i = 0; i < vect.size(); ++i) {
				norm_sq += vect[i] * vect[i];
			}
			return norm_sq;
		}
	};
}

int main() {
	int n = 10;
	vector<double> iterate(10, 1.0);
	int iteration_count = gradient_descent(iterate, Objective<>(), Retraction<>());
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
