CFLAGS=-std=gnu99 -Wall -O0

all: main main-secuencial

main: main.c
	$(CC) $(CCFLAGS) -fopenmp -o main main.c

main-secuencial: main.c
	$(CC) $(CCFLAGS) -o main-secuencial main.c


clean:
	rm -rf main
