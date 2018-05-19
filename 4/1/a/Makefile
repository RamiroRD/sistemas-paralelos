CC=mpicc
CFLAGS=-std=c99 -O0 -Wall -pedantic-errors -ggdb

all: p2p

imprimir_rank: 

run: p2p
	mpirun -np 4 p2p


clean:
	rm -f p2p
