LIBLAPACKPATH = /project/gavoth/ynhan/local/lapack-3.5.0/
LAPACKEPATH = /project/gavoth/ynhan/local/lapack-3.5.0/lapacke/include
OPT = -O3 -Wall

CFLAGS = $(OPT) -I$(LAPACKEPATH)
LDFLAGS = $(OPT) -L$(LIBLAPACKPATH)
LIBS = -llapacke -llapack -lblas -lgfortran
CC = g++

test.x : test.o
	$(CC) $(LDFLAGS) -o $@ test.o $(LIBS)

test.o : test.cpp
	$(CC) $(CFLAGS) -c test.cpp
