ROOTLIBS   = $(shell root-config --libs)
ROOTCFLAGS = $(shell root-config --cflags)

CXX        = gcc
CXXFLAGS   = $(ROOTCFLAGS) -Wall -g
LD         = gcc
LDFLAGS    = -Wall -g
LIBS       = $(ROOTLIBS) -lEG

L1RatesAndPrescales: L1RatesAndPrescales.cpp L1Tree.h L1Tree.cpp
	$(CXX) $(CXXFLAGS) $(LIBS) L1Tree.cpp L1RatesAndPrescales.cpp -o L1RatesAndPrescales

HltRates: HltRates.cpp HltTree.h HltTree.cpp
	$(CXX) $(CXXFLAGS) $(LIBS) HltTree.cpp HltRates.cpp -o HltRates


all: L1RatesAndPrescales HltRates
