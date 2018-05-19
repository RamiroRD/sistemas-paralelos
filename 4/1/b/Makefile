CC=mpicc
CFLAGS=-std=c99 -O0 -Wall -pedantic-errors -ggdb

all: colectivas

imprimir_rank: 

run: colectivas
	mpirun -np 4 colectivas


clean:
	rm -f colectivas
