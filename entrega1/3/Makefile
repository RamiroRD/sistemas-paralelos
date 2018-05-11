CFLAGS=-std=gnu99 -Wall -O0

all: pares pares-secuencial

pares: pares.c
	$(CC) $(CFLAGS) -fopenmp -o pares pares.c

pares-secuencial: pares.c
	$(CC) $(CFLAGS) -o pares-secuencial pares.c

clean:
	rm -rf pares pares-secuencial

