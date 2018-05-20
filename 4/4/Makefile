CC=mpicc
CFLAGS=-std=c99 -O0 -Wall -pedantic-errors -ggdb

all: secuencial paralelo

secuencial: secuencial.c
	gcc $(CFLAGS) -o secuencial secuencial.c

paralelo: 

clean:
	rm -f secuencial paralelo
