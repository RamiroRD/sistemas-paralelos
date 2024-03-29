#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int t;
uint64_t n;
int *V;

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

int main(int argc, char **argv)
{
	/* Tiempos */
	double ti, tf;

	if (argc != 3) {
		printf("pares T N\n");
		return -1;
	}

	/* Cantidad de hilos */
	t = atoi(argv[1]);
	/* Dimensión del arreglo */
	n = strtoll(argv[2], NULL, 10);
	if (errno == ERANGE) {
		printf("Fuera de rango.\n");
		return -1;
	}


	fprintf(stderr, "Alocando %f MiB...\n", sizeof(int) * n / 1048576.0);
	V = malloc(sizeof(int) * n);
	if (!V) {
		perror("malloc");
		return -1;
	}


	for (uint64_t i = 0; i < n; i++)
		V[i] = i;

#ifdef _OPENMP
	omp_set_num_threads(t);
#endif

	ti = dwalltime();

	uint64_t pares = 0;
#pragma omp parallel for reduction(+:pares)
	for (uint64_t i = 0; i < n; i++)
		if((V[i] & 1) == 0)
			pares++;

	tf = dwalltime();

	printf("T = %f [s]\n", tf - ti);
	fprintf(stderr, "Cantidad de pares: %lu\n", pares);

	free(V);

	return 0;
}

