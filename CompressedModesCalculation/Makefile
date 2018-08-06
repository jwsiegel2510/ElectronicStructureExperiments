#This makefile builds the LaplaceCompressedModes program.
#
## Usage:
#
## make LaplaceCompressedModes 
#

CC = g++
OPT = -lm -std=c++11

LaplaceCompressedModes : compressed_modes_objective.h laplace_compressed_modes.cpp
	$(CC) -o LaplaceCompressedModes.exe laplace_compressed_modes.cpp $(OPT)
