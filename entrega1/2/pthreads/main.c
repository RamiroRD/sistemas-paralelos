#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>

int t, n;
double *A, *B, *C, *D, *E, *F, *U, *L;
double *AT, *BT, *CT, *ET, *FT, *UT;
/* Intermedios */
double *AA, *AAC;
double *LB, *LBE;
double *DU, *DUF;
/* Sumas de U, L y B junto a sus mutexes*/ 
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

#warning No implementado
/* 
 * Traspone src, dejando el resultado en dst. src y dst deben ser distintos.
 */
void transpose(double * restrict dst, const double * restrict src, int n, int t, int id)
{
}

void multiply(double *C, const double * restrict B, const double * restrict A, int n, int t, int id)
{
	int slice = n / t;
	for (int i = id * slice; i < slice * (id + 1); i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[j * n + k];
}

#warning No implementado
void transpose_upper(double * restrict dst, const double * restrict src, int n, int t, int id)
{
}

#warning No implementado
/*
 * Suma A y B, dejando el resultado en C. A y B deben estar almacenadas por
 * filas. C puede ser A o B.
 */
void add(double *C, const double *B, const double *A, int n, int t, int id)
{
}

#warning No implementado
/* 
 * Multiplica cada elemento de A por factor. Deja el resultado directamente en
 * A
 */
void scale(double *A, double factor, int n, int t, int id)
{
}

#warning No implementado
/*
 * Suma todos los elementos de A. Añade la suma a la variable res atómicamente,
 * usando el mútex mutex.
 */
void sum(const double *A, double *res, pthread_mutex_t *mutex, int n, int id)
{
}

#warning No implementado
/*
 * Calcula el producto A * B donde A es una matriz triangular inferior.
 * El resultado se almacena en C. C != B != A.
 */
void multiply_ll(double * restrict C, const double * restrict A, const double * restrict B, int n, int t, int id)
{
}

#warning No implementado
/*
 * Calcula el producto A * B donde B es una matriz triangular superior.
 * El resultado se almacena en C. C != B != A.
 */
void multiply_ru(double * restrict C, const double * restrict A, const double * restrict B, int n, int t, int id)
{
}

void *worker(void *idp)
{
	int id = *(int*) idp;

	/* Promedio de u */
	sum(U, &up, &up_mutex, n, id);
	pthread_barrier_wait(&barrier);
	/* Promedio de l */
	sum(L, &lp, &lp_mutex, n, id);
	pthread_barrier_wait(&barrier);
	/* Promedio de b */
	sum(L, &bp, &lp_mutex, n, id);
	pthread_barrier_wait(&barrier);

	/* Transpuestas */
	transpose(AT, A, n, t, id);
	transpose(BT, B, n, t, id);
	transpose(CT, C, n, t, id);
	transpose(ET, E, n, t, id);
	transpose(FT, F, n, t, id);
	transpose_upper(UT, U, n, t, id);
	pthread_barrier_wait(&barrier);

	/*  AA */
	AA = C; /* Reutilizamos el espacio de C */
	multiply(C, A, AT, n, t, id);
	pthread_barrier_wait(&barrier);

	/* AAC */
	AAC = A; /* Reutilizamos A */
	multiply(AAC, AA, C, n, t, id);
	
	/* ulAAC */
	scale(C, (up + lp) / (n * n * n * n), n * n, t, id);
	pthread_barrier_wait(&barrier);

	/* LB */
	LB = C;  /* Reutilizamos el espacio de C de vuelta */
	multiply_ll(LB, L, B, n, t, id);
	pthread_barrier_wait(&barrier);

	/* LBE */
	LBE = B; /* Reutilizamos el espacio de B */
	multiply(LBE, LB, E, n, t, id);
	pthread_barrier_wait(&barrier);

	/* DU */
	DU = C; /* Reutilizamos el espacio de C */
	multiply_ru(DU, D, U, n, t, id);
	pthread_barrier_wait(&barrier);

	/* DUF */
	DUF = E; /* Reutilizamos el espacio de C */
	multiply(DUF, DU, F, n, t, id);
	pthread_barrier_wait(&barrier);

	/* b/(LBE + DUF) en el espacio de C */
	add(C, LBE, DUF, n, t, id);
	scale(C, bp / (n * n) , n, t, id);

	/* Resultado final en A */
	add(A, AAC, C, n, t, id);


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


	int ids[t];
	pthread_t threads[t];

	for (int i = 0; i < t; i++)
		ids[i] = i;

	/*
	 * Inicializamos los mútex.
	 */
	pthread_mutex_init(&up_mutex, NULL);
	pthread_mutex_init(&lp_mutex, NULL);
	pthread_mutex_init(&bp_mutex, NULL);

	/*
	 * Reservamos memoria para las matrices A, B, C, D, E, F, L, U y para
	 * la representación ordenada por columnas de A, B, C, E, F, U.
	 */
	A = malloc(sizeof(double) * n * n);
	B = malloc(sizeof(double) * n * n);
	C = malloc(sizeof(double) * n * n);
	D = malloc(sizeof(double) * n * n);
	E = malloc(sizeof(double) * n * n);
	F = malloc(sizeof(double) * n * n);
	L = malloc(sizeof(double) * (n * (n + 1) / 2));
	U = malloc(sizeof(double) * (n * (n + 1) / 2));
	AT = malloc(sizeof(double) * n * n);
	BT = malloc(sizeof(double) * n * n);
	CT = malloc(sizeof(double) * n * n);
	ET = malloc(sizeof(double) * n * n);
	FT = malloc(sizeof(double) * n * n);
	UT = malloc(sizeof(double) * (n * (n + 1) / 2));

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
	pthread_mutex_destroy(&up_mutex);
	pthread_mutex_destroy(&lp_mutex);
	pthread_mutex_destroy(&bp_mutex);

	return 0;
}

