CC=mpicc
CFLAGS=-std=c99 -Wall -pedantic-errors

all: imprimir_rank

imprimir_rank: 

run: imprimir_rank
	mpirun -np 4 imprimir_rank


clean:
	rm -f imprimir_rank
