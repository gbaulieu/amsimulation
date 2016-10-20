UNAME := $(shell uname)

SRC=src

FLAG=-O3 -Wall -std=c++11 -Werror=type-limits
INC =-I `root-config --incdir` -I ${SRC}
LIBS =-L ${ROOTSYS}/lib

module_positions_generator: module_positions_generator.o ModuleData.o
	g++ -o module_positions_generator ModuleData.o module_positions_generator.o ${LIBS} -lCore -lCint -lRIO -lHist -lTree -lMatrix -lGpad

ModuleData.o: ModuleData.h ModuleData.cc
	g++ -c $(LIBS) ModuleData.cc

clean:	
	rm -rf ./*.o;rm -f ./*~;rm -f ./*#

module_positions_generator.o:./module_positions_generator.cc
	g++ -c ${FLAG} ${INC} ./module_positions_generator.cc

