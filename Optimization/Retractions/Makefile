
# This makefile builds and runs tests for the libraries in this directory.
#
# Usage:
#
# make test --builds and runs the test programs
#

CC = g++
OPT = -lm -std=c++11
INCLUDES = -I./ -I../../EigenLib
test : .stiefel_cayley_retraction_test

.stiefel_cayley_retraction_test : stiefel_cayley_retraction.h stiefel_cayley_retraction_test.cpp
	$(CC) -o stiefel_cayley_retraction_test.exe stiefel_cayley_retraction_test.cpp $(OPT) $(INCLUDES)
	./stiefel_cayley_retraction_test.exe
