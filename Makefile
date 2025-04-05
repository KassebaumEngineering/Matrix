# Makefile for test of matrix class
#

CXX=g++
CXXFLAGS=-g -O -Wall -Wextra -std=c++11 # Using CXX, CXXFLAGS, adding warnings and C++11 standard
LIBS=-lm # Removed -lg++, g++ handles stdc++ automatically

all: simple twolevel

twolevel: twolevel.o Matrix.o
	$(CXX) -o twolevel $(CXXFLAGS) twolevel.o Matrix.o $(LIBS)

twolevel.o: twolevel.C Matrix.H
	$(CXX) -c $(CXXFLAGS) twolevel.C

simple: simple.o Matrix.o
	$(CXX) -o simple $(CXXFLAGS) simple.o Matrix.o $(LIBS)

simple.o: simple.C Matrix.H
	$(CXX) -c $(CXXFLAGS) simple.C

matrixtest: matrixtest.o Matrix.o TimeUse.o TimeValue.o
	$(CXX) -o matrixtest $(CXXFLAGS) matrixtest.o Matrix.o TimeUse.o TimeValue.o $(LIBS)

matrixtest.o: matrixtest.C Matrix.H TimeUse.H
	$(CXX) -c $(CXXFLAGS) matrixtest.C

Matrix.o: Matrix.C
	$(CXX) -c -O2 Matrix.C # Keeping -O2 for Matrix.o for now

TimeUse.o: TimeUse.C TimeUse.H TimeValue.H
	$(CXX) -c $(CXXFLAGS) TimeUse.C

TimeValue.o: TimeValue.C TimeValue.H
	$(CXX) -c $(CXXFLAGS) TimeValue.C

clean:
	rm -f *.o simple matrixtest twolevel

