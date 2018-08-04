/* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains a templated implementation of preconditioned gradient descent on a manifold.
 *  Passed as template parameters are a class representing points on the manifold, 
 *  a class representing dual tangent vectors (i.e. gradients), a class evaluating 
 *  the gradient and objective, a class performing the retraction (i.e. moving 
 *  on the manifold), and a class which applies the preconditioner. The class 
 *  which performs the retraction must also expose a method which evaluates the 
 *  norm of a dual tangent vector and method which applies the preconditioner must
 *  also return the norm of the preconditioned gradient.
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

#ifndef _OPTIMIZATION_METHODS_PRECONDITIONED_GRADIENT_DESCENT__
#define _OPTIMIZATION_METHODS_PRECONDITIONED_GRADIENT_DESCENT__

#include<cmath>

namespace optimization {
namespace methods {

template<class P, class V, template<class, class> class Objective, template<class, class> class Retraction, template<class, class> class Preconditioner>
int preconditioned_gradient_descent(P& iterate, Objective<P,V> objective, Retraction<P,V> retraction, Preconditioner<P,V> preconditioner, double gtol = 1e-4) {
	const static double rho = 5e-1; // Parameter for the Armijo condition.
	const static double eta = .5; // Parameter determining the factor by which the step decreases if Armijo condition isn't met.
	double step_size = 1e-1; // Initial step size.

	P temporary_iterate(iterate);
	V grad;
	double obj = objective.evaluate(iterate);

	int it = 0;
	while(1) {
		// Evaluate gradient.
		objective.evaluate_grad(grad, iterate);

		// Calculate the squared norm of the gradient and check stopping condition.
		double grad_norm_sq = retraction.norm_sq(grad, temporary_iterate);
		if (sqrt(grad_norm_sq) < gtol) {
			return it;
		}

		// Apply preconditioner
		double preconditioned_grad_norm_sq = preconditioner.precondition(grad, iterate);

		// Calculate the retraction. Use a line search to determine the step size.
		temporary_iterate = iterate;
		bool done_increase = false;
		while (!done_increase) {
			step_size /= eta;
			retraction.retract(iterate, grad, -1.0 * step_size);
			double F = objective.evaluate(iterate);
			if (F <= obj - rho * step_size * preconditioned_grad_norm_sq) {
				step_size /= eta;
				iterate = temporary_iterate;
			} else {
				done_increase = true;
			}
		}
		bool done_decrease = false;
		while (!done_decrease) {
			retraction.retract(iterate, grad, -1.0 * step_size);
			double F = objective.evaluate(iterate);
                        if (F <= obj - rho * step_size * preconditioned_grad_norm_sq) {
                                done_decrease = true;
				obj = F;
                        } else {
                                step_size *= eta;
				iterate = temporary_iterate;
                        }
		}
		++it;
	}
}

} // namespace methods
} // namespace optimization

#endif
