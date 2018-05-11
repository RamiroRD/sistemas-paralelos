CFLAGS= -std=gnu99 -Wall -O0

all: producto producto-secuencial

producto: producto.c
	$(CC) $(CFLAGS) -fopenmp -o producto producto.c

producto-secuencial: producto.c
	$(CC) $(CFLAGS) -o producto-secuencial producto.c

clean:
	rm -rf producto producto-secuencial

