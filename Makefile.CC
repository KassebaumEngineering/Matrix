#
# Makefile for test of Matrix class
#
#  Revision:
#  $Id: Makefile.CC,v 1.3 1993/11/27 00:20:18 jak Exp $
#
#  History:
#  $Log: Makefile.CC,v $
#  Revision 1.3  1993/11/27 00:20:18  jak
#  Matrix Class has been ported for use with the AT&T cfront compiler version 3
#  (with templates).   -jak
#
# Revision 1.2  1993/11/15  20:29:36  jak
# Corrections and fixes.  Works now with GCC2.5.3 and Libg++2.5.1 -jak
#
# =====================================

CC=CC
CFLAGS=-O5
LIBS=-lMatrixCC -lC -lm -liostream
INSTALL_LIB=/usr/local/lib
INSTALL_INC=/usr/local/include

SRC=Makefile.gcc Makefile.CC matrixtest.C Matrix.C Matrix.H \
Linked_List_Template.H TimeValue.C TimeValue.H TimeUse.C TimeUse.H

HEADERS=Matrix.H Linked_List_Template.H TimeValue.H TimeUse.H
OBJ= Matrix.o TimeValue.o TimeUse.o
EXES=matrixtest
OTHER= matrixtest.o libMatrixCC.a

all: libMatrixCC.a matrixtest 

matrixtest: matrixtest.o libMatrixCC.a
	$(CC) -o matrixtest $(CFLAGS) matrixtest.o -L. $(LIBS)

matrixtest.o: matrixtest.C Matrix.H Linked_List_Template.H
	$(CC) -c $(CFLAGS) matrixtest.C

Matrix.o: Matrix.C Matrix.H Linked_List_Template.H
	$(CC) -c $(CFLAGS) Matrix.C

TimeValue.o: TimeValue.C TimeValue.H
	$(CC) -c $(CFLAGS) TimeValue.C

TimeUse.o: TimeUse.C TimeUse.H Linked_List_Template.H TimeValue.H
	$(CC) -c $(CFLAGS) TimeUse.C

libMatrixCC.a: $(OBJ)
	-rm -f libMatrixCC.a
	ar cr libMatrixCC.a $(OBJ)
	ranlib libMatrixCC.a

install:
	rm -f $(INSTALL_LIB)/libMatrixCC.a
	cp libMatrixCC.a $(INSTALL_LIB)
	ranlib $(INSTALL_LIB)/libMatrixCC.a
	for i in $(HEADERS) ; do \
	  cp $$i $(INSTALL_INC)/$$i ; \
	done

clean:
	rm -f $(OBJ) $(EXES) $(OTHER)

dist: $(SRC)
	mkdir Matrix
	cp -p $(SRC) Matrix
	for i in $(SRC) ; do \
	  cp $$i Matrix/$$i ; \
	done
	tar crf - Matrix | compress -c > Matrix.tar.Z
	rm -rf Matrix
	uuencode Matrix.tar.Z Matrix.tar.Z > Matrix.uu
	touch dist
