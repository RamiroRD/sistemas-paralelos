CFLAGS=-std=gnu99 -Wall -O0

all: main main-secuencial

main: main.c
	$(CC) $(CFLAGS) -fopenmp -o main main.c

main-secuencial: main.c
	$(CC) $(CFLAGS) -o main-secuencial main.c


clean:
	rm -rf main main-secuencial
