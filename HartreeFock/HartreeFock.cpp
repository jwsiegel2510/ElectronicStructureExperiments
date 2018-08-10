/* This file is part of an open source library for electronic structure calculations.
 *
 *  It contains code which performs a Hartree-Fock calculation based on parameters passed 
 *  via an xml file. A file containing the orbital integrals for a desired basis set must be 
 *  provided as well.
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
#include "RealOrbitalIntegralsData.h"
#include "Eigen/Dense"
#include "accelerated_gradient_descent.h"
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
	using ::optimization::methods::accelerated_gradient_descent;
	using ::optimization::retractions::StiefelCayleyRetraction;
	using ::optimization::preconditioners::StiefelQuadraticPreconditioner;
	using ::SCC_Q::RealOrbitalIntegralsData;

	template <class P = MatrixXd, class V = MatrixXd> class HartreeFockObjective {
	private:
		const RealOrbitalIntegralsData& integrals_data;
		int basis_count;
	
	public:
		HartreeFockObjective() = delete;
		HartreeFockObjective(const RealOrbitalIntegralsData& input_data) : integrals_data(input_data) {
			basis_count = integrals_data.getOrbitalBasisCount();
		}

		// The input is expected to have 2*basis_count rows (due to spin effects).
		double evaluate(const& MatrixXd iterate) {
			double output = 0;
			// Evaluate the single particle energy.
			for (int j = 0; j < iterate.cols(); ++j) {
				for (int i = 0; i < basis_count; ++i) {
					for (int k = 0; k < basis_count; ++k) {
						output += integrals_data(i,k) * (iterate(2*i,j) * iterate(2*k,j) + iterate(2*i + 1,j) * iterate(2*k + 1,j);
					}
				}
			}
			// Evaluate the electron-electron interaction energy.
			for (int j1 = 0; j1 < iterate.cols(); ++j1) {
				for (int j2 = j1 + 1; j2 < iterate.cols(); ++j2) {
					for (int i1 = 0; i1 < basis_count; ++i1) {
						for (int i2 = 0; i2 < basis_count; ++i2) {
							for (int k1 = 0; k1 < basis_count; ++k1) {
								for (int k2 = 0; k2 < basis_count; ++k2) {
									output += integrals_data(i1,i2,k1,k2) * 
											(iterate(2*i1,j1) * iterate(2*i2,j1) + iterate(2*i1+1,j1) * iterate(2*i2+1,j1)) * 
											(iterate(2*k1,j2) * iterate(2*k2,j2) + iterate(2*k1+1,j2) * iterate(2*k2+1,j2));
									output -= integrals_data(i1,k1,i2,k2) * 
											(iterate(2*i1,j1) * iterate(2*k1,j2) + iterate(2*i1+1,j1) * iterate(2*k1+1,j2)) * 
											(iterate(2*i2,j1) * iterate(2*k2,j2) + iterate(2*i2+1,j1) * iterate(2*k2+1,j2));
								}
							}
						}
					}
				}
			}
			return output;
		}

		void evaluate_grad(MatrixXd& grad, const MatrixXd& iterate) {
			grad = MatrixXd::Zero(input.rows(), input.cols());
                        for (int j = 0; j < iterate.cols(); ++j) {
                                for (int i = 0; i < basis_count; ++i) {
                                        for (int k = 0; k < basis_count; ++k) {
                                                output += integrals_data(i,k) * (iterate(2*i,j) * iterate(2*k,j) + iterate(2*i + 1,j) * iterate(2*k + 1,j);
						grad(2*i,j) += integrals_data(i,k) * iterate(2*k,j);
						grad(2*i+1,j) += integrals_data(i,k) * iterate(2*k+1,j);
                                                grad(2*k,j) += integrals_data(i,k) * iterate(2*i,j);
                                                grad(2*k+1,j) += integrals_data(i,k) * iterate(2*i+1,j);
                                        }
                                }
                        }
                        for (int j1 = 0; j1 < iterate.cols(); ++j1) {
                                for (int j2 = j1 + 1; j2 < iterate.cols(); ++j2) {
                                        for (int i1 = 0; i1 < basis_count; ++i1) {
                                                for (int i2 = 0; i2 < basis_count; ++i2) {
                                                        for (int k1 = 0; k1 < basis_count; ++k1) {
                                                                for (int k2 = 0; k2 < basis_count; ++k2) {
									grad(2*i1,j1) += integrals_data(i1,i2,k1,k2) * iterate(2*i2,j1) * (iterate(2*k1,j2) * iterate(2*k2,j2) + iterate(2*k1+1,j2) * iterate(2*k2+1,j2));
									grad(2*i2,j1) += integrals_data(i1,i2,k1,k2) * iterate(2*i1,j1) * (iterate(2*k1,j2) * iterate(2*k2,j2) + iterate(2*k1+1,j2) * iterate(2*k2+1,j2));
									grad(2*i1+1,j1) += integrals_data(i1,i2,k1,k2) * iterate(2*i2+1,j1)) * (iterate(2*k1,j2) * iterate(2*k2,j2) + iterate(2*k1+1,j2) * iterate(2*k2+1,j2));
									grad(2*i2+1,j1) += integrals_data(i1,i2,k1,k2) * iterate(2*i1+1,j1) * (iterate(2*k1,j2) * iterate(2*k2,j2) + iterate(2*k1+1,j2) * iterate(2*k2+1,j2));
									grad(2*k1,j2) += integrals_data(i1,i2,k1,k2) * iterate(2*k2,j2) * (iterate(2*i1,j1) * iterate(2*i2,j1) + iterate(2*i1+1,j1) * iterate(2*i2+1,j1));
									grad(2*k2,j2) += integrals_data(i1,i2,k1,k2) * iterate(2*k1,j2) * (iterate(2*i1,j1) * iterate(2*i2,j1) + iterate(2*i1+1,j1) * iterate(2*i2+1,j1));
									grad(2*k1+1,j2) += integrals_data(i1,i2,k1,k2) * iterate(2*k2+1,j2) * (iterate(2*i1,j1) * iterate(2*i2,j1) + iterate(2*i1+1,j1) * iterate(2*i2+1,j1));
									grad(2*k2+1,j2) += integrals_data(i1,i2,k1,k2) * iterate(2*k1+1,j2) * (iterate(2*i1,j1) * iterate(2*i2,j1) + iterate(2*i1+1,j1) * iterate(2*i2+1,j1));
                                                                        output -= integrals_data(i1,k1,i2,k2) *
                                                                                        (iterate(2*i1,j1) * iterate(2*k1,j2) + iterate(2*i1+1,j1) * iterate(2*k1+1,j2)) *
                                                                                        (iterate(2*i2,j1) * iterate(2*k2,j2) + iterate(2*i2+1,j1) * iterate(2*k2+1,j2));
									grad(2*i1,j1) -= integrals_data(i1,k1,i2,k2) * iterate(2*k1,j2) * (iterate(2*i2,j1) * iterate(2*k2,j2) + iterate(2*i2+1,j1) * iterate(2*k2+1,j2));
									grad(2*k1,j2) -= integrals_data(i1,k1,i2,k2) * iterate(2*i1,j1) * (iterate(2*i2,j1) * iterate(2*k2,j2) + iterate(2*i2+1,j1) * iterate(2*k2+1,j2));
									grad(2*i1+1,j1) -= integrals_data(i1,k1,i2,k2) * iterate(2*k1+1,j2) * (iterate(2*i2,j1) * iterate(2*k2,j2) + iterate(2*i2+1,j1) * iterate(2*k2+1,j2));
									grad(2*k1+1,j2) -= integrals_data(i1,k1,i2,k2) * iterate(2*i1+1,j1) * (iterate(2*i2,j1) * iterate(2*k2,j2) + iterate(2*i2+1,j1) * iterate(2*k2+1,j2));
									grad(2*i2,j1) -= integrals_data(i1,k1,i2,k2) * iterate(2*k2,j2) * (iterate(2*i1,j1) * iterate(2*k1,j2) + iterate(2*i1+1,j1) * iterate(2*k1+1,j2));
									grad(2*k2,j2) -= integrals_data(i1,k1,i2,k2) * iterate(2*i2,j1) * (iterate(2*i1,j1) * iterate(2*k1,j2) + iterate(2*i1+1,j1) * iterate(2*k1+1,j2));
									grad(2*i2+1,j1) -= integrals_data(i1,k1,i2,k2) * iterate(2*k2+1,j2) * (iterate(2*i1,j1) * iterate(2*k1,j2) + iterate(2*i1+1,j1) * iterate(2*k1+1,j2));
									grad(2*k2+1,j2) -= integrals_data(i1,k1,i2,k2) * iterate(2*i2+1,j1) * (iterate(2*i1,j1) * iterate(2*k1,j2) + iterate(2*i1+1,j1) * iterate(2*k1+1,j2));
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
		}
	};
}

int main() {
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
	double tol = paramList.getParameterValue("tolerance", "Parameters");
	int k = paramList.getParameterValue("electron_count", "Parameters");
	string integral_data_file = paramList.getParameterValue("integral_data", "Parameters");

	// Load Orbital Integral data.
	ReadOrbitalIntegralsData integrals_data;
	integrals_data.inputOrbitalIntegralData(integral_data_file);
	int n = integrals_data.getOrbitalBasisCount();

	// Perform Calculation.
	MatrixXd iterate(2*n,k);
	HartreeFockObjective<> objective(integrals_data);
	StiefelCayleyRetraction<MatrixXd, MatrixXd> retraction;
	retraction.generate_random_point(iterate);
	std::cout << "Iteration Count: " << accelerated_gradient_descent(iterate, objective, retraction, tol) << "\n";
	
	// Output ground state energy.
	std::cout objective.evaluate(iterate);
}
