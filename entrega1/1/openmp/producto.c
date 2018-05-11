#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int t, n;
double *A, *B, *C;

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
		printf("producto T n\n");
		return -1;
	}

	/* Cantidad de hilos */
	t = atoi(argv[1]);
	/* Dimensión de bloque (en elementos) */
	n = atoi(argv[2]);

	A = malloc(sizeof(double) * n * n);
	B = malloc(sizeof(double) * n * n);
	C = malloc(sizeof(double) * n * n);

	/* Inicializamos C en cero */
	memset(C, 0, sizeof(double) * n * n);

	for (int i = 0; i < n * n; i++)
		A[i] = 1.0;

#ifdef _OPENMP
	omp_set_num_threads(t);
#endif

	ti = dwalltime();

	/*
	 * El operando de mano derecha tiene que ser A pero almacenado por
	 * columnas.
	 */
#pragma omp parallel for default(none) shared(n, B, A)
	for (int i = 0; i < n; i++)	
		for (int j = 0; j < n; j++)	
			B[i * n + j] = A[j * n + i];

	/*
	 * Hacemos la multiplicación convencional.
	 */
#pragma omp parallel for default(none) shared(n, C, B, A)
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[j * n + k];
	tf = dwalltime();

	bool ok = true;
	for (int i = 0; i < n * n; i++)
		ok = ok && C[i] == n;

	fprintf(stderr, ok ? "OK\n" : "ERROR\n");
	printf("T = %f [s]\n", tf - ti);

	free(A);
	free(B);
	free(C);

	return 0;
}

