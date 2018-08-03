# This makefile builds and runs tests for the libraries in this directory.
#
# Usage:
#
# make test --builds and runs the test programs
#

CC = g++
OPT = -lm -std=c++11
test : .gradient_descent_test .accelerated_gradient_descent_test .preconditioned_gradient_descent_test .preconditioned_accelerated_gradient_descent_test

.gradient_descent_test : gradient_descent.h test_lib.h gradient_descent_test.cpp
	$(CC) -o gradient_descent_test.exe gradient_descent_test.cpp $(OPT)
	./gradient_descent_test.exe

.accelerated_gradient_descent_test : accelerated_gradient_descent.h test_lib.h accelerated_gradient_descent_test.cpp
	$(CC) -o accelerated_gradient_descent_test.exe accelerated_gradient_descent_test.cpp $(OPT)
	./accelerated_gradient_descent_test.exe

.preconditioned_gradient_descent_test : preconditioned_gradient_descent.h test_lib.h preconditioned_gradient_descent_test.cpp
	$(CC) -o preconditioned_gradient_descent_test.exe preconditioned_gradient_descent_test.cpp $(OPT)
	./preconditioned_gradient_descent_test.exe
	 

.preconditioned_accelerated_gradient_descent_test : preconditioned_accelerated_gradient_descent.h test_lib.h preconditioned_accelerated_gradient_descent_test.cpp
	$(CC) -o preconditioned_accelerated_gradient_descent_test.exe preconditioned_accelerated_gradient_descent_test.cpp $(OPT)
	./preconditioned_accelerated_gradient_descent_test.exe
