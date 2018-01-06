#
# Makefile for stepps2 STAN simulations.
#

#COMP: INTEL, GNU
#COMP = INTEL
COMP ?= GNU

STAN     ?= $(abspath $(PREDICTION)/../cmdstan-2.17.1)
STANC    ?= $(STAN)/bin/stanc
STANMAIN ?= $(STAN)/src/cmdstan/main.cpp
STANLIBS ?= $(STAN)/stan/lib/stan_math/lib

ifeq ($(COMP), GNU)
  CXX = g++
endif

ifeq ($(COMP), INTEL)
  CXX = icl
endif

# include local config
include $(wildcard *.defs)

CXXFLAGS += -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -DBOOST_DATE_TIME_NO_LIB
CXXFLAGS += -I$(STAN)/stan/src
CXXFLAGS += -I$(STAN)/src
CXXFLAGS += -I$(STAN)/stan/lib/stan_math
CXXFLAGS += -I$(STANLIBS)/eigen_3.3.3
CXXFLAGS += -I$(STANLIBS)/boost_1.64.0
CXXFLAGS += -I$(STANLIBS)/cvodes_2.9.0/include

#
# rules
#

all:

%.cpp: %.stan
	$(STANC) --o=$*.cpp $*.stan

%.exe: %.cpp
	cat $*.cpp $(STANMAIN) > $*_tmp.cpp
	$(CXX) -o $*.exe $(CXXFLAGS) $*_tmp.cpp $(LDFLAGS)

clean:
	rm -rf build *.exe

.PHONY: clean
.PRECIOUS: %.cpp
