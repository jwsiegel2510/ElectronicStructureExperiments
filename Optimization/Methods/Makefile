# This makefile builds and runs tests for the libraries in this directory.
#
# Usage:
#
# make test --builds and runs the test programs
#

CC = g++
OPT = -lm -std=c++11
test : .gradient_descent_test

.gradient_descent_test: gradient_descent_test.exe
	./gradient_descent_test.exe

gradient_descent_test.exe : gradient_descent.h gradient_descent_test.cpp
	$(CC) -o gradient_descent_test.exe gradient_descent_test.cpp $(OPT)
	 
