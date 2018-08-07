# This makefile builds the LaplaceCompressedModes program.
#
# Jonathan W. Siegel 2018
#
# To build the executable using this script use the
# command
#
# make release 
#
# or 
#
# make debug
#
# To see a verbose output of the build add VERBOSE=1 to the make invocation commands above 

SHELL=/bin/sh

# Specifying the location of the MakeScripts *.mk files

SCRIPTS_PATH=../MakeScripts

# Use BaseCommonConfig.mk if it exists otherwise use BaseCommonConfig_Default.mk 

ifneq ("$(wildcard $(SCRIPTS_PATH)/BaseCommonConfig.mk)","")
	include $(SCRIPTS_PATH)/BaseCommonConfig.mk
else
	include $(SCRIPTS_PATH)/BaseCommonConfig_Default.mk
endif

# Specify that level 4 optimization is to be used.

CXXFLAGS+= -O4

CPPfiles+= ./laplace_compressed_modes.cpp

INCLUDES+= -I./ -I../Optimization/Methods -I../Optimization/Retractions 
INCLUDES+= -I../XML-ParameterList -I../EigenLib -I../SCC-Utility

# Library specification 

LIBS     += -lXML_ParameterList  
LIB_PATH += -L../XML-ParameterList/lib

# Specification of the release and debug directories and executables.

RELEASE_DIR= ./_releaseLaplaceCompressedModes

DEBUG_DIR= ./_debugLaplaceCompressedModes

RELEASE_EXEC= LaplaceCompressedModes.exe

DEBUG_EXEC=LaplaceCompressedModesDebug.exe

include $(SCRIPTS_PATH)/ExecutableMake.mk
