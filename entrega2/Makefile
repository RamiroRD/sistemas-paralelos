CC=cc
CFLAGS=-std=c99 -O0 -Wall -pedantic-errors -ggdb -Wno-unknown-pragmas

all: secuencial mpi hibrido

mpi: paralelo.c
	mpicc $(CFLAGS) -o mpi paralelo.c

hibrido: paralelo.c
	mpicc $(CFLAGS) -fopenmp -o hibrido paralelo.c

clean:
	rm -f secuencial mpi hibrido
