# -*- Makefile -*- 

#
# This Makefile is intended for compiling cc under delphes
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#Tao Huang, email: taohuang@physics.tamu.edu
#include Makefile.arch
#LDLIBS=`root-config --glibs`
#ROOTLIBS='-lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats -lHistFactory'
ROOTINCLUDE=$(shell root-config --glibs 2> /dev/null) 
#Other_INCLUDE_if_necessary= 
LDFLAGS=-Wl,--no-as-needed
LDFLAGS+=`root-config --ldflags`
INCLUDE = $(ROOTINCLUDE)
#
# C++ flags
# 
CXX=g++
CXXFLAGS=-pedantic -ansi -Wall -Wno-long-long -Wno-format -Werror=uninitialized -Werror=delete-non-virtual-dtor -O2
CXXFLAGS+=`root-config --cflags`
#src file
Example_SRC = testMake.cpp 
Example_EXE = testMake.exe

#DiHiggsAnalyzer
TestHME_SRC = testHME.cc  heavyMassEstimator.h heavyMassEstimator.cc
TestHME_EXE = testHME.exe
# test print out
#var=TEST
#PHONY:
#	@echo $(INCLUDE)
	#echo $(LIB)
#default : DiHiggs
all: testHME
testHME:
	$(CXX) $(CXXFLAGS) -g $(TestHME_SRC) $(LDFLAGS) $(INCLUDE) -o $(TestHME_EXE)	

clean:
	rm -f *.exe
