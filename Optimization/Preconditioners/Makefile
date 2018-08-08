# This makefile builds and runs tests for the libraries in this directory.
#
# Usage:
#
# make test --builds and runs the test programs
#

CC = g++
OPT = -lm -std=c++11
INCLUDES = -I./ -I../../EigenLib -I../Retractions
test : .stiefel_quadratic_preconditioner_test

.stiefel_quadratic_preconditioner_test : stiefel_quadratic_preconditioner.h stiefel_quadratic_preconditioner_test.cpp ../Retractions/stiefel_cayley_retraction.h
	$(CC) -o stiefel_quadratic_preconditioner_test.exe stiefel_quadratic_preconditioner_test.cpp $(OPT) $(INCLUDES)
	./stiefel_quadratic_preconditioner_test.exe

