#
# Makefile for STEPPS calibration model
# 

#COMP: INTEL, GNU
#COMP = INTEL
COMP = GNU

#STAN     ?= $(STEPPS2)/opt/stan210
STAN     ?= $(HOME)/Documents/paleon/stepps2/opt/stan210
STANC    ?= $(STAN)/bin/stanc

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

CXXFLAGS += -I$(STAN)/src
CXXFLAGS += -I$(STAN)/lib/eigen_3.2.0
CXXFLAGS += -I$(STAN)/lib/boost_1.54.0
CXXFLAGS += -I$(STAN)/lib/gtest_1.7.0
#CXXFLAGS += -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -pipe -lpthread

#
# rules
#

all:

%.cpp: %.stan
	$(STANC) --o=$*.cpp $*.stan
#	$(STEPPS2)/addtimer.pl build/$*.cpp

%.exe: %.cpp
	$(CXX) -o $*.exe $(CXXFLAGS) $*.cpp $(LDFLAGS)

clean:
	rm -rf build *.exe

.PHONY: clean
.PRECIOUS: %.cpp
