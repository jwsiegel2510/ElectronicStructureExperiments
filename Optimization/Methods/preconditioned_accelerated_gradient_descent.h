/* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains a templated implementation of preconditioned accelerated gradient 
 *  descent on a manifold.
 *
 *  Passed as template parameters are a class representing points on the manifold, 
 *  a class representing dual tangent vectors (i.e. gradients), a class evaluating 
 *  the gradient and objective, a class performing the retraction (i.e. moving 
 *  on the manifold), and a class which applies the preconditioner. 
 *  The class which performs the retraction must also expose a 
 *  method which evaluates the norm of a dual tangent vector and a method which 
 *  extrapolates on the manifold. The method which applies the preconditioner must also
 *  return the norm of the preconditioned vector.
 *
 *  The templated classes which implement the objective, retraction, and preconditioner must take
 *  the point and dual tangent vector classes as template arguments. They may also
 *  take additional template arguments as long as the point and dual tangent vector
 *  classes are the only template arguments shared between them.
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

#ifndef _OPTIMIZATION_METHODS_GRADIENT_DESCENT__
#define _OPTIMIZATION_METHODS_GRADIENT_DESCENT__

#include<cmath>

namespace optimization {
namespace methods {
template<class P, class V, template<class, class, typename ... > class Objective, template<class, class, typename ... > class Retraction, 
template<class, class, typename ... > class Preconditioner, typename ... AdditionalArgsObj, typename ... AdditionalArgsRtr, typename ... AdditionalArgsPC>
int preconditioned_accelerated_gradient_descent(P& iterate, Objective<P, V, AdditionalArgsObj ... > objective, 
							    Retraction<P, V, AdditionalArgsRtr ... > retraction, 
							    Preconditioner<P, V, AdditionalArgsPC ... > preconditioner, double gtol = 1e-4) {	
	const static double restart_rho = 1e-2; // Parameter determining the sufficient decrease which triggers a restart.
	const static double rho = 5e-1; // Parameter for the Armijo condition (must be set to at least .5 so that the step size is leq 1/L).
	const static double eta = .5; // Determines the factor by which the step size decreases if the Armijo condition isn't met.
	double step_size = 1e-1; // Initial step size.

	P y_iterate(iterate);
	P temporary_iterate(y_iterate);
	V grad;

	double obj;
	int it = 0;
	int k = 0;
	while (1) {
		// Take a gradient step starting from the y_iterate.
		obj = objective.evaluate(y_iterate);
		objective.evaluate_grad(grad, y_iterate);

		double grad_norm_sq = retraction.norm_sq(grad, y_iterate);
		if (sqrt(grad_norm_sq) < gtol) {
			iterate = y_iterate; return it;
		}

		// Apply Preconditioner
		double preconditioned_grad_norm_sq = preconditioner.precondition(grad, iterate);

		// Decrease step size until the Armijo condition is met (allow an increase if k = 0).
		temporary_iterate = y_iterate;
		bool done_increase = (k > 0);
		while (!done_increase) {
			retraction.retract(y_iterate, grad, -1.0 * step_size);
			double F = objective.evaluate(y_iterate);
			if (F <= obj - rho * step_size * preconditioned_grad_norm_sq) {
				step_size /= eta;
				y_iterate = temporary_iterate;
			} else {
				done_increase = true;
			}
		}
		bool done_decrease = false;
		while (!done_decrease) {
			retraction.retract(y_iterate, grad, -1.0 * step_size);
			double F = objective.evaluate(y_iterate);
                        if (F <= obj - rho * step_size * preconditioned_grad_norm_sq) {
                                done_decrease = true;
				obj = F;
                        } else {
                                step_size *= eta;
				y_iterate = temporary_iterate;
                        }
		}

		// Restart momentum if there is not a sufficient decrease in the objective.
		if (objective.evaluate(y_iterate) > objective.evaluate(iterate) - restart_rho * step_size * grad_norm_sq) {
			y_iterate = iterate;
			k = 0;
		} else {
			temporary_iterate = y_iterate;
			retraction.extrapolate(y_iterate, iterate, temporary_iterate, k/(k + 3.0));
			iterate = temporary_iterate;
			++k;
		}
		++it;
	}
}

} // namespace methods
} // namespace optimization
	
#endif
