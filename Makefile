# Makefile for test of matrix class
#

CC=/usr/local/bin/gcc
CFLAGS=-g -O
LIBS=-lg++ -lm

all: simple twolevel

twolevel: twolevel.o Matrix.o
	$(CC) -o twolevel $(CFLAGS) twolevel.o Matrix.o $(LIBS)

twolevel.o: twolevel.C Matrix.H
	$(CC) -c $(CFLAGS) twolevel.C

simple: simple.o Matrix.o
	$(CC) -o simple $(CFLAGS) simple.o Matrix.o $(LIBS)

simple.o: simple.C Matrix.H
	$(CC) -c $(CFLAGS) simple.C

matrixtest: matrixtest.o Matrix.o
	$(CC) -o matrixtest $(CFLAGS) matrixtest.o Matrix.o $(LIBS)

matrixtest.o: matrixtest.C Matrix.H
	$(CC) -c $(CFLAGS) matrixtest.C

Matrix.o: Matrix.C
	$(CC) -c -O2 Matrix.C

clean:
	rm -f *.o simple matrixtest twolevel

