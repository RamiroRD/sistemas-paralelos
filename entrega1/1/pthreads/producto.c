#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>
#include <pthread.h>

int t, n;
double *A, *B, *C;
pthread_barrier_t barrier;

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

void *worker(void *idp)
{
	int id = *(int*) idp;
	int slice = n / t;

	pthread_barrier_wait(&barrier);
	for (int i = id * slice; i < slice * (id + 1); i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[j * n + k];
	pthread_barrier_wait(&barrier);
	return NULL;
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

	int ids[t];
	pthread_t threads[t];

	for (int i = 0; i < t; i++)
		ids[i] = i;

	/* Inicializamos A */
	for (int i = 0; i < n * n; i++)
		A[i] = 1.0;

	/* Inicializamos C en cero */
	memset(C, 0, sizeof(double) * n * n);

	/*
	 * El operando de mano derecha tiene que ser A pero almacenado por
	 * columnas.
	 */
	for (int i = 0; i < n; i++)	
		for (int j = 0; j < n; j++)	
			B[i * n + j] = A[j * n + i];

	/*
	 * Creamos la barrera y los hilos.
	 * No se utiliza el primer elemento a propósito para dejar el código
	 * más legible.
	 */
	pthread_barrier_init(&barrier, NULL, t);
	for (int i = 1; i < t; i++)
		pthread_create(threads + i, NULL, worker, ids + i);

	/*
	 * Hacemos la multiplicación convencional.
	 */
	ti = dwalltime();
	worker(ids + 0);
	tf = dwalltime();

	bool ok = true;
	for (int i = 0; i < n * n; i++)
		ok = ok && C[i] == n;

	fprintf(stderr, ok ? "OK\n" : "ERROR\n");
	printf("T = %f [s]\n", tf - ti);

	free(A);
	free(B);
	free(C);
	pthread_barrier_destroy(&barrier);
	
	return 0;
}

