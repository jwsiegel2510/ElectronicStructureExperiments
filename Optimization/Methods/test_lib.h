 /* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains simple quadratic functions for use in the tests of the various optimization methods. 
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

#ifndef _OPTIMIZATION_METHODS_TEST_LIB__
#define _OPTIMIZATION_METHODS_TEST_LIB__

#include<vector>

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

	void extrapolate(P& output, const P& start, const P& end, double overshoot) {
		output.resize(start.size());
		for (int i = 0; i < start.size(); ++i) {
			output[i] = start[i] + (1 + overshoot) * (end[i] - start[i]);
		}
	}
};


#endif
