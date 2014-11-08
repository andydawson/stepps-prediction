#
# Makefile for stepps2 STAN simulations.
#

#COMP: INTEL, GNU
#COMP = INTEL
COMP = GNU

STAN     ?= $(STEPPS2)/opt/stan240
STANC    ?= $(STAN)/bin/stanc
MAIN     ?= $(STAN)/src/cmdstan/main.cpp

ifeq ($(COMP), GNU)
  CXX      = g++
 # CXXFLAGS = -O3 -fopenmp 
  LDFLAGS  ?= -L$(STAN)/bin -lstan -lrt	
endif

ifeq ($(COMP), INTEL)
  CXX      = icpc
  CXXFLAGS = -O3 -openmp #-ipo -no-prec-div -openmp 
  LDFLAGS  ?= -L$(STAN)/bin -lstan -lrt
endif

# include local config
include $(wildcard *.defs)

CXXFLAGS += -I$(STAN)/stan/src
CXXFLAGS += -I$(STAN)/stan/lib/eigen_3.2.0
CXXFLAGS += -I$(STAN)/stan/lib/boost_1.54.0
CXXFLAGS += -I$(STAN)/stan/lib/gtest_1.7.0
#CXXFLAGS += -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -pipe -lpthread

#
# rules
#

all:

%.cpp: %.stan
	$(STANC) --o=$*.cpp $*.stan

%.exe: %.cpp
	$(CXX) -o $*.exe $(CXXFLAGS) $(MAIN) -include $*.cpp $(LDFLAGS)

clean:
	rm -rf build *.exe

.PHONY: clean
.PRECIOUS: %.cpp
