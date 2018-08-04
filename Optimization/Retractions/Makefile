
# This makefile builds and runs tests for the libraries in this directory.
#
# Usage:
#
# make test --builds and runs the test programs
#

CC = g++
OPT = -lm -std=c++11
test : .stiefel_cayley_retraction_test

.stiefel_cayley_retraction_test : stiefel_cayley_retraction.h stiefel_cayley_retraction_test.cpp
	$(CC) -o stiefel_cayley_retraction_test.exe stiefel_cayley_retraction_test.cpp $(OPT)
	./stiefel_cayley_retraction_test.exe
