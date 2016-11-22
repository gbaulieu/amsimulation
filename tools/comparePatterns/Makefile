UNAME := $(shell uname)

SRC=src

CMSSW_ENABLED=true

ifeq ($(CMSSW_ENABLED),true)
	FLAG=-O3 -Wall -std=c++11 -Werror=type-limits
	INC =-I `root-config --incdir` -I `scram tool tag boost INCLUDE` -I ${SRC}
	LIBS =-L ${ROOTSYS}/lib -L `scram tool tag boost LIBDIR`
else
	FLAG=-O3 -Wall -std=c++11 -Werror=type-limits
	INC =-I `root-config --incdir` -I ${SRC}
	LIBS =-L ${ROOTSYS}/lib
endif
BOOSTLIBS = -lboost_serialization -lboost_program_options -lboost_iostreams

comparePatterns: comparePatterns.o Event.o
	g++ -o comparePatterns Event.o comparePatterns.o ${LIBS} ${BOOSTLIBS} -lCore -lCint -lRIO -lHist -lTree -lMatrix -lGpad

clean:	
	rm -rf ./*.o;rm -f ./*~;rm -f ./*#

comparePatterns.o:./comparePatterns.cc
	g++ -c ${FLAG} ${INC} ./comparePatterns.cc

Event.o: Event.cc Event.h
	g++ -c ${FLAG} ${INC} ./Event.cc