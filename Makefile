#kmeans : c++ implementation of the kmeans algorithm (http://en.wikipedia.org/wiki/K-means_clustering)
#for finding clusters of stable microstates when studying ligand migration in proteins
# 
#Copyright (c) 2015, Florent Hédin, Pierre-André Cazade, and the University of Basel.
#All rights reserved.
# 
#The 3-clause BSD license is applied to this software.
#See LICENSE.txt

#################################################################
########################   MakeVars   ###########################
#################################################################

CXX=g++
# CXX=icpc

CXX_OPT= -std=c++11 -I "./include" -march=native -Wall -Wextra -pedantic -O3
# CXX_OPT= -std=c++11 -I "./include" -march=native -Wall -Wextra -pedantic -O0 -g

LD_LIB=

LD_OPT=-lm -lblas -llapack

MKDIR=mkdir -p ./obj

TARGET=kmeans

SRC=$(wildcard ./src/*.cpp)

OBJ=$(patsubst ./src/%.cpp,./obj/%.o,$(SRC))

#################################################################
########################   Makefile   ###########################
#################################################################

all:$(TARGET)
	@echo "Compilation Success"

$(TARGET):Makefile

./obj/%.o:./src/%.cpp
	@$(MKDIR)
	$(CXX) $(CXX_OPT) -c $< -o $@

$(TARGET):$(OBJ)
	$(CXX) $(CXX_OPT) $(LD_LIB) $(OBJ) -o $@ $(LD_OPT)

clean:
	rm -f $(TARGET) ./obj/*.o
