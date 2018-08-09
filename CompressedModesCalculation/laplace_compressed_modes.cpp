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

#include "XML_ParameterListArray.h"
#include "Eigen/Dense"
#include "preconditioned_accelerated_gradient_descent.h"
#include "stiefel_cayley_retraction.h"
#include "stiefel_quadratic_preconditioner.h"
#include "CmdOptionUtility.h"
#include "compressed_modes_objective.h"
#include <cstdio>
#include <iostream>
#include <fstream>

namespace {
	using ::Eigen::MatrixXd;
	using ::Eigen::IOFormat;
	using ::optimization::methods::preconditioned_accelerated_gradient_descent;
	using ::optimization::retractions::StiefelCayleyRetraction;
	using ::optimization::preconditioners::StiefelQuadraticPreconditioner;
	using ::compressed_modes::CompressedModesObjective;

	class LaplaceOperator { // 1D Laplacian operator on Eigen matrices.
	private:
		int n;

	public:

		LaplaceOperator(int n_) : n(n_) {}	

		MatrixXd apply(const MatrixXd& input) const {
			MatrixXd output(input);
			for (int j = 0; j < input.cols(); ++j) {
				for (int i = 0; i < input.rows(); ++i) {
					if (i == 0) {
						output(i,j) = (2 * input(i,j) - input(i+1,j)) * n;
					} else if (i == (input.rows() - 1)) {
						output(i,j) = (input(i,j) - input(i-1,j)) * n;
					} else {
						output(i,j) = (2 * input(i,j) - input(i+1,j) - input(i-1,j)) * n;
					}
				}
			}
			return output;
		}

		void invert(MatrixXd& input) const {
			if (input.rows() <= 1) {
				input /= (2 * n); return;
			}
			for (int j = 0; j < input.cols(); ++j) {
				for (int i = input.rows() - 2; i >= 0; --i) {
					input(i,j) += input(i+1,j);
				}
				for (int i = 1; i < input.rows(); ++i) {
					input(i,j) += input(i-1,j);
				}
				for (int i = 0; i < input.rows(); ++i) input(i,j) /= n;
			}
		}
	};
}

int main(int argc, char** argv) {
	// Get input XML filename and output filename.
	string parameterFileName;
	string outputFileName;
	CmdOptionUtility optionUtility;
	parameterFileName = optionUtility.getCmdOption(argc,argv,"-f");
	outputFileName = optionUtility.getCmdOption(argc,argv,"-o");
	if (parameterFileName.empty() || outputFileName.empty()) {
		printf("XML filename must be passed with -f flag and output filename must be passed with -o flag.\n");
		return -1;
	}

	// Initialize parameters from XML file.
	XML_ParameterListArray paramList;
	paramList.setAbortOnErrorFlag();
	paramList.initialize(parameterFileName.c_str());
	double mu = paramList.getParameterValue("mu", "Parameters");
	double epsilon = paramList.getParameterValue("epsilon", "Parameters");
	double tol = paramList.getParameterValue("tolerance", "Parameters");
	int n = paramList.getParameterValue("GridPoints", "Parameters");
	int k = paramList.getParameterValue("ModesCount", "Parameters");

	// Perform Calculation.
	MatrixXd iterate(n,k);
	LaplaceOperator op_(n);
	CompressedModesObjective<MatrixXd, MatrixXd, LaplaceOperator> objective(op_, mu / n, epsilon);
	StiefelQuadraticPreconditioner<MatrixXd, MatrixXd, LaplaceOperator> preconditioner(op_);
	StiefelCayleyRetraction<MatrixXd, MatrixXd> retraction;
	retraction.generate_random_point(iterate);
	std::cout << "Iteration Count: " << preconditioned_accelerated_gradient_descent(iterate, objective, retraction, preconditioner, tol) << "\n";

	// Output Result.
	const IOFormat fmt(-1, 1, "\t", " \n ", "(", ")", "\n", "\n");
	std::ofstream out;
	out.open(outputFileName);
	out << iterate.format(fmt) << "\n";
}
