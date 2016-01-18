UNAME := $(shell uname)

SRC=src

CMSSW_ENABLED=false

ifeq ($(CMSSW_ENABLED),true)
	FLAG=-O3 -Wall -std=c++11 -Werror=type-limits
	INC =-I `root-config --incdir` -I `scram tool tag boost INCLUDE` -I ${SRC}
	LIBS =-L ${ROOTSYS}/lib -L `scram tool tag boost LIBDIR`
else
	FLAG=-O3 -Wall -Werror=type-limits
	INC =-I `root-config --incdir` -I ${SRC}
	LIBS =-L ${ROOTSYS}/lib
endif
BOOSTLIBS = -lboost_serialization -lboost_program_options -lboost_iostreams

mergeSectors: mergeSectors.o
	g++ -o mergeSectors mergeSectors.o ${LIBS} ${BOOSTLIBS} -lCore -lCint -lRIO -lHist -lTree -lMatrix -lGpad

clean:	
	rm -rf ./*.o;rm -f ./*~;rm -f ./*#

mergeSectors.o:./mergeSectors.cc
	g++ -c ${FLAG} ${INC} ./mergeSectors.cc

