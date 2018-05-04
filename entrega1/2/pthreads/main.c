#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>

int t, n;
double *A, *B, *C, *D, *E, *F, *U, *L;
/* Promedios de U, L y B junto a sus mutexes*/ 
double up, lp, bp;
pthread_mutex_t up_mutex, lp_mutex, bp_mutex;
pthread_barrier_t barrier;

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

#warning Implementar
void transpose(double * restrict dst, const double * restrict src, int n, int t, int id)
{
}

#warning Implementar
void multiply(double *C, const double *B, const double *A, int n, int t, int id)
{
	int slice = n / t;
	for (int i = id * slice; i < slice * (id + 1); i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[j * n + k];
}

#warning Implementar
void add(double *C, const double *B, const double *A, int n, int id)
{
}

#warning Implementar
void scale(double *A, double factor, int n, int id)
{
}

#warning Implementar
void sum(const double *A, double *res, pthread_mutex_t *mutex,int n, int id)
{
}

void *worker(void *idp)
{
	int id = *(int*) idp;

	/* Promedio de u */
	pthread_barrier_wait(&barrier);
	/* Promedio de l */
	pthread_barrier_wait(&barrier);
	/* Promedio de b */
	pthread_barrier_wait(&barrier);

	/* Transposición de A */
	pthread_barrier_wait(&barrier);
	/* Transposición de B */
	pthread_barrier_wait(&barrier);
	/* Transposición de C */
	pthread_barrier_wait(&barrier);
	/* Transposición de E */
	pthread_barrier_wait(&barrier);
	/* Transposición de F */
	pthread_barrier_wait(&barrier);
	/* Transposición de U */
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

	/*
	 * Inicializamos la barrera y arrancamos los hilos.
	 */
	pthread_barrier_init(&barrier, NULL, t);
	ti = dwalltime();
	for (int i = 0; i < t; i++)
		pthread_create(threads + i, NULL, worker, ids + i);
	for (int i = 0; i < t; i++)
		pthread_join(threads[i], NULL);
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

