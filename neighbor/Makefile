LIBLAPACKPATH = /project/gavoth/ynhan/local/lapack-3.5.0/
LAPACKEPATH = /project/gavoth/ynhan/local/lapack-3.5.0/lapacke/include
CBLASPATH = /project/gavoth/ynhan/local/CBLAS/include
LIBCBLAS = /project/gavoth/ynhan/local/CBLAS/lib/
OPT = -O3 -Wall

CFLAGS = -I$(LAPACKEPATH) -I$(CBLASPATH)
LDFLAGS = -L$(LIBLAPACKPATH) -L$(LIBCBLAS)
LIBS = -lblas -llapacke -llapack -lcblas -lgfortran
CC = g++ $(OPT)

test.x : test.o engine.o fg_atoms.o cg_sites.o matrix_C.o matrix_M.o matrix_N.o mapping.o neighbor.o
	$(CC) -o test.x *.o $(LDFLAGS) $(LIBS)

test.o : test.cpp
	$(CC) -c test.cpp

mapping.o : mapping.cpp
	$(CC) -c mapping.cpp

engine.o : engine.cpp
	$(CC) $(CFLAGS) -c engine.cpp

fg_atoms.o : fg_atoms.cpp
	$(CC) -c fg_atoms.cpp

cg_sites.o : cg_sites.cpp
	$(CC) -c cg_sites.cpp

matrix_C.o : matrix_C.cpp
	$(CC) -c matrix_C.cpp

matrix_M.o : matrix_M.cpp
	$(CC) -c matrix_M.cpp

matrix_N.o : matrix_N.cpp
	$(CC) -c matrix_N.cpp

neighbor.o : neighbor.cpp
	$(CC) -c neighbor.cpp

clean :
	rm *.o
