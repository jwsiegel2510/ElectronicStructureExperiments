/* This file is part of an open source library for electronic structure calculations.
 *  
 *  It contains a templated implementation of a retraction on the Stiefel manifold
 *  based on the Cayley transform (alternatively on the (1,1) Pade approximation of
 *  of the exponential. The classes passed as parameters are assumed to represent
 *  points on the manifold and dual tangent vectors. 
 *  
 *  As such, instances of these classes must essentitally be matrices, 
 *  supporting linear algebra operations, including solving linear systems.
 *  The synthax is assumed to be the same as for Eigen matrices.
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

#ifndef _OPTIMIZATION_RETRACTIONS_STIEFEL_CAYLEY_RETRACTION__
#define _OPTIMIZATION_RETRACTIONS_STIEFEL_CAYLEY_RETRACTION__

#include <cstdlib>
#include <ctime>
#include <cmath>

namespace optimization {
namespace retractions {

template<class P, class V> class StiefelCayleyRetraction {
private:
	P U;
	P VM;
	V X;
public:
        void retract(P& iterate, const V& direction, double dt);
        void extrapolate(P& output, const P& point_one, const P& point_two, double beta);
        double norm_sq(const V& grad, const P& iterate);
	void generate_random_point(P& iterate);
};

template<class P, class V> void StiefelCayleyRetraction<P,V>::retract(P& iterate, const V& direction, double dt) {
	U.resize(iterate.rows(), 2 * iterate.cols());
	VM.resize(iterate.rows(), 2 * iterate.cols());
	U.block(0, 0, iterate.rows(), iterate.cols()) = -0.5 * dt * direction;
	U.block(0, iterate.cols(), iterate.rows(), iterate.cols()) = 0.5 * dt * iterate;
	VM.block(0, 0, iterate.rows(), iterate.cols()) = iterate;
	VM.block(0, iterate.cols(), iterate.rows(), iterate.cols()) = direction;
	iterate = iterate - 2.0 * U * (P::Identity(2 * iterate.cols(), 2 * iterate.cols()) + VM.transpose() * U).partialPivLu().solve(VM.transpose() * iterate);
}

template<class P, class V> double StiefelCayleyRetraction<P,V>::norm_sq(const V& grad, const P& iterate) {
        X.resize(grad.rows(), grad.cols());
	X = grad;
	X -= iterate * (grad.transpose() * iterate);
        return (X.transpose() * X).trace();
}

template<class P, class V> void StiefelCayleyRetraction<P,V>::extrapolate(P& output, const P& point_one, const P& point_two, double beta) {
	X = 2.0 * (P::Identity(point_one.cols(), point_one.cols()) + point_two.transpose() * point_one).partialPivLu().solve(point_two.transpose()).transpose();
	output = point_one;
	retract(output, X, 1.0 + beta);
}

template<class P, class V> void StiefelCayleyRetraction<P,V>::generate_random_point(P& iterate) {
	srand(time(NULL));
	for (int j = 0; j < iterate.cols(); ++j) {
		for (int i = 0; i < iterate.rows(); ++i) {
			double theta = 2 * M_PI * ((double) rand()) / (RAND_MAX);
			double r = ((double) rand()) / (RAND_MAX);
			iterate(i,j) = sqrt(-log(r)) * sin(theta);
		}
	}
	iterate = iterate.householderQr().householderQ() * P::Identity(iterate.rows(), iterate.cols());
}

} // namespace retractions
} // namespace optimization

#endif
