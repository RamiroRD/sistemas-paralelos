#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>

/* Defines para acceso a matrices triangulares */
#define U_FIL(i,j) (i * n + j - (i * (i + 1)) / 2)
#define U_COL(i,j) (i + (j * (j + 1)) / 2)
#define L_FIL(i,j) (j + (i * (i + 1)) / 2)

/* Cantidad de hilos y dimensión de matrices */
int t, n;

/* Matrices entrada y alguna de sus transpuestas */
double *A, *B, *C, *D, *E, *F, *U, *L;
double *AT, *BT, *CT, *ET, *FT, *UT;

/* Matrices intermedias. No se reservan memoria exclusiva para estos*/
double *AA, *AAC;
double *LB, *LBE;
double *DU, *DUF;

/* Sumas de U, L y B */
double up, lp, bp;


double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

/*
 * Traspone src, dejando el resultado en dst. src y dst deben ser distintos.
 */
void transpose(double * restrict dst, const double * restrict src, int n)
{
  #pragma omp for
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			dst[i * n + j] = src[j * n + i];
}

void multiply(double *C, const double * restrict B, const double * restrict A, int n)
{
	/* Multiplicación convencional fila * columna */
  #pragma omp for
  for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
      double partial = 0;
			for (int k = 0; k < n; k++)
				partial += A[i * n + k] * B[j * n + k];
			C[i * n + j] = partial;
}
#TODO hacerlo
void transpose_upper(double * restrict dst, const double * restrict src, int n)
{
	#pragma omp for
	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++)
			dst[U_COL(i,j)] = src[U_FIL(i,j)];

}

/*
 * Suma A y B, dejando el resultado en C. A y B deben estar almacenadas por
 * filas. C puede ser A o B.
 */
void add(double *C, const double *A, const double *B, int n)
{
  #pragma omp for
	for (int i = 0; i < n * n; i++)
		C[i] = A[i] + B[i];
}

/*
 * Multiplica cada elemento de A por factor. Deja el resultado directamente en
 * A.
 */
void scale(double *A, double factor, int dim)
{
  #pragma omp for
	for (int i = 0; i < dim; i++)
		A[i] *= factor;
}

/*
 * Suma todos los elementos de A. Añade la suma a la variable res atómicamente,
 * usando el mútex mutex.
 */
void sum(const double *A, double *res, pthread_mutex_t *mutex, int n)
{
	double partial = 0;
  #pragma omp for reduction(+:partial)
	for (int i = 0; i < n * n; i++)
		partial += A[i];
	*res = partial;
}

/*
 * Calcula el producto A * B donde A es una matriz triangular inferior ordenada
 * por filas.
 * El resultado se almacena en C. C != B != A.
 */
void multiply_ll(double * restrict C, const double * restrict A, const double * restrict B, int n)
{
  /* Multiplicación convencional fila * columna */
  #pragma omp for
  for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
    {
			double partial = 0;
			/* Cuando k > i los elementos de A valen 0, por lo tanto se deja de sumar. */
			for (int k = 0; k <= i; k++)
			 partial += A[L_FIL(i,k)] * B[j * n + k];
			C[i * n + j] = partial;
    }

}

/*
 * Calcula el producto A * B donde B es una matriz triangular superior ordenada
 * por columnas.
 * El resultado se almacena en C. C != B != A.
 */
void multiply_ru(double * restrict C, const double * restrict A, const double * restrict B, int n)
{
	/* Multiplicación convencional fila * columna */
  #pragma omp for
  for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
    {
      double partial = 0;
			/* Cuando k > j los elementos de B valen 0, por lo tanto se deja de sumar. */
			for (int k = 0; k <= j; k++)
				 partial += A[i * n + k] * B[U_COL(k,j)];
			C[i * n + j] = partial;
    }

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
	 * Inicializamos arrancamos los hilos.
	 */

  omp_set_num_threads(t);

	ti = dwalltime();

  #pragma omp parallel {

  /* Promedio de u */
	sum(U, &up, n);
	#pragma omp barrier
	/* Promedio de l */
	sum(L, &lp, n);
	#pragma omp barrier
	/* Promedio de b */
	sum(B, &bp, n);
	#pragma omp barrier

	/* Transpuestas */
	transpose(AT, A, n);
	transpose(BT, B, n);
	transpose(CT, C, n);
	transpose(ET, E, n);
	transpose(FT, F, n);
	transpose_upper(UT, U, n);
	#pragma omp barrier

	/*  AA */
	AA = C; /* Reutilizamos el espacio de C */
	multiply(C, A, AT, n);
	#pragma omp barrier

	/* AAC */
	AAC = A; /* Reutilizamos A */
	multiply(AAC, AA, CT, n);

	/* ulAAC */
	scale(AAC, (up + lp) / (n * n), n * n);
	#pragma omp barrier

	/* LB */
	LB = C;  /* Reutilizamos el espacio de C de vuelta */
	multiply_ll(LB, L, BT, n);
	#pragma omp barrier

	/* LBE */
	LBE = B; /* Reutilizamos el espacio de B */
	multiply(LBE, LB, ET, n);
	#pragma omp barrier

	/* DU */
	DU = C; /* Reutilizamos el espacio de C */
	multiply_ru(DU, D, UT, n);
	#pragma omp barrier

	/* DUF */
	DUF = E; /* Reutilizamos el espacio de E*/
	multiply(DUF, DU, FT, n);
	#pragma omp barrier

	/* b/(LBE + DUF) en el espacio de C */
	add(C, LBE, DUF, n);
	scale(C, bp / (n * n) , n);

	/* Resultado final en A */
	add(A, AAC, C, n);

  }
	tf = dwalltime();

	printf("T = %f [s]\n", tf - ti);

	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	free(L);
	free(U);
	free(AT);
	free(BT);
	free(CT);
	free(ET);
	free(FT);
	free(UT);

	return 0;
}
