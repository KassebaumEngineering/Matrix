#
# Makefile for test of Matrix class
#
#  Revision:
#  $Id: Makefile.CC,v 1.2 1993/11/15 20:29:36 jak Exp $
#
#  History:
#  $Log: Makefile.CC,v $
#  Revision 1.2  1993/11/15 20:29:36  jak
#  Corrections and fixes.  Works now with GCC2.5.3 and Libg++2.5.1 -jak
#*/
# =====================================

CC=CC
CFLAGS=-O -DATT_CC
LIBS=-lC -lm

SRC=Makefile.gcc Makefile.CC matrixtest.C Matrix.C Matrix.H simple.C twolevel.C

all: matrixtest simple twolevel

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
	$(CC) -c $(CFLAGS) Matrix.C

clean:
	rm -f *.o simple matrixtest twolevel

dist: $(SRC)
	mkdir Matrix
	cp Makefile.* Matrix
	cp *.C Matrix
	cp *.H Matrix
	tar crf - Matrix | compress -c > Matrix.tar.Z
	rm -rf Matrix
	uuencode Matrix.tar.Z Matrix.tar.Z > Matrix.uu
	rm -rf Matrix.tar.Z
	touch dist

