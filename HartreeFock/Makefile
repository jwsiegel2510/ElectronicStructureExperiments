
#
# This is the "meta" makefile that builds the executable and it's dependents
# by invoking their respective makescripts.
#
# Usage :
#
# make release                   === builds release versions of all executables
# make debug                     === builds debug versions of all executables
#
# make clean                     === removes object files and executables
# make cleanall                  === removes object files, executables and temporary directories
#
# Jan. 18, 2018
#
SHELL=/bin/sh

# Location of Makescripts files

MAKESCRIPTS_Dir=../MakeScripts

# Use BaseCommonConfig.mk if it exists, otherwise use BaseCommonConfig_Default.mk 

ifneq ("$(wildcard $(MAKESCRIPTS_Dir)/BaseCommonConfig.mk)","")
	include $(MAKESCRIPTS_Dir)/BaseCommonConfig.mk
else
	include $(MAKESCRIPTS_Dir)/BaseCommonConfig_Default.mk
endif


ifeq ($(MAKECMDGOALS),release)
BUILD_TYPE=release
endif
ifeq ($(MAKECMDGOALS),debug)
BUILD_TYPE=debug
endif

release : .Starting   .XML_ParameterListLib  .HartreeFock   .Finished 

debug   : .Starting   .XML_ParameterListLib_debug .HartreeFock_debug .Finished  

.Starting :
	###########################################################
	# Build Started
	########################################################### 
	$(QUIET) echo "Compilation Date : " `date` 
ifeq ($(OpenMP),0)
	$(QUIET) echo "Executable Type  :  Single-threaded" 
else
	$(QUIET) echo "Executable Type  :  OpenMP based multi-threaded" 
endif

 
.XML_ParameterListLib :
	###########################################################
	# .XML_ParameterListLib
	###########################################################
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk --directory=../XML-ParameterList release 
.XML_ParameterListLib_debug :
	###########################################################
	# .XML_ParameterListLib_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk --directory=../XML-ParameterList debug


.HartreeFock :
	###########################################################
	# .OrbitalIntegralsMain_AU
	###########################################################
	$(QUIET)$(MAKE) -f HartreeFock.mk  release 
.HartreeFock_debug :
	###########################################################
	# .OrbitalIntegralsMain_AU_debug 
	########################################################### 
	$(QUIET)$(MAKE) -f HartreeFock.mk  debug
	


.Finished :
	###########################################################
	# Build Completed
	########################################################### 
	
.PHONY : cleanall clean
cleanall  : clean
	rm -rf ./_*
clean: 
	$(QUIET)$(MAKE) -f XML_ParameterListLib.mk --directory=../XML-ParameterList   clean 
	$(QUIET)$(MAKE) -f HarteeFock.mk clean
