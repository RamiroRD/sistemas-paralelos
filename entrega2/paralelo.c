#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <unistd.h>
#include <sys/time.h>
#include <fcntl.h>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define TAG 1
#define INLINE __attribute__((always_inline)) inline

/* Defines para acceso a matrices triangulares */
#define U_FIL(i,j) ((i) * (n) + (j) - ((i) * ((i) + 1)) / 2)
#define U_COL(i,j) ((i) + ((j) * ((j) + 1)) / 2)
#define L_FIL(i,j) ((j) + ((i) * ((i) + 1)) / 2)

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;

	return sec;
}

/*
 * Suma A y B, dejando el resultado en C. A y B deben estar almacenadas por
 * filas. C puede ser A o B.
 */
INLINE
void add(double *C, const double *A, const double *B, int n, int t, int rank)
{
	const int slice = n * n / t;

#pragma omp parallel for
	for (int i = 0; i < slice; i++)
		C[i] = A[i] + B[i];
}

/*
 * Multiplica cada elemento de A por factor. Deja el resultado directamente en
 * A.
 */
INLINE
void scale(double *A, double factor, int n, int t, int rank)
{
	const int slice = n * n / t;

#pragma omp parallel for
	for (int i = 0; i < slice; i++)
		A[i] *= factor;
}

/*
 * Suma todos los elementos de A. Añade la suma a la variable res atómicamente,
 * usando el mútex mutex.
 */
INLINE
double sum(const double *A, int n, int t, int rank)
{
	const int slice = n * n / t;

	double s = 0;

#pragma omp parallel for reduction(+:s)
	for (int i = 0; i < slice; i++)
		s += A[i];

	return s;
}

INLINE
void multiply(double *C, const double *restrict A,
	      const double *restrict B, int n, int t, int rank)
{
	const int slice = n / t;
	const int offset = rank * slice;
	double c;

	/* Multiplicación convencional fila * columna */
#pragma omp parallel for private(c)
	for (int i = offset; i < (rank + 1) * slice; i++) {
		for (int j = 0; j < n; j++) {
			c = 0;
			for (int k = 0; k < n; k++) {
				c += A[(i - offset) * n + k] * B[j * n + k];
			}
			C[(i - offset) * n + j] = c;
		}
	}

}

/*
 * Calcula el producto A * B donde A es una matriz triangular inferior ordenada
 * por filas.
 * A la vez, calcula la suma de los elementos de A.
 * El resultado se almacena en C. C != B != A.
 */
INLINE
double multiply_ll(double *restrict C, const double *restrict A,
		 const double *restrict B, int n, int t, int rank)
{
	const int slice = n / t;
	const int offset = rank * slice;
	double c;
	double s = 0;

#pragma omp parallel for schedule(dynamic, 16) private(c,s)
	for (int i = offset; i < (rank + 1) * slice; i++) {
		for (int j = 0; j < n; j++) {
			c = 0;
			for (int k = 0; k <= i; k++)
				c += A[L_FIL(i - offset, k)] * B[j * n + k];

			C[(i - offset) * n + j] = c;
			if (i >= j)
				s += A[L_FIL(i - offset,j)];
		}
	}
	return s;
}

/*
 * Calcula el producto A * B donde B es una matriz triangular superior ordenada
 * por columnas.
 * A la vez, calcula la suma de los elementos de B.
 * El resultado se almacena en C. C != B != A.
 */
INLINE
double multiply_ru(double *restrict C, const double *restrict A,
		 const double *restrict B, int n, int t, int rank)
{
	const int slice = n / t;
	const int offset = rank * slice;
	double c;
	double s = 0;

#pragma omp parallel for schedule(dynamic, 16) private(c,s)
	for (int i = offset; i < (rank + 1) * slice; i++) {
		for (int j = 0; j < n; j++) {
			c = 0;
			for (int k = 0; k <= j; k++)
				c += A[(i - offset) * n + k] * B[U_COL(k, j)];

			C[(i - offset) * n + j] = c;
			if (i <= j)
				s += B[U_COL(i - offset, j)];
		}
	}

	return s;
}

void wait_for_gdb()
{
	int x = 0;
	while (!x)
		sleep(1);
}

/*
 * Compara las matrices A y B con una tolerancia epsilon. A y B deben estar
 * almacenadas de la misma forma.
 */
void error(double *A, double *B, int n, double *resultados)
{
	double sum = 0;
	double res = 0;
	for (int i = 0; i < n * n; i++) {
		double d = (A[i] - B[i]);
		d *= d;
		res += d;
		sum += B[i];
	}
	resultados[0] = res / (n * n);
	resultados[1] = sum / (n * n);
}

void load(const char *filename, double *A, double *B, double *C, double *D,
		double *L, double *U, double *res, int n)
{
	int fd = open(filename, O_RDONLY);

	if (fd < 1) {
		perror("load");
		exit(-1);
	}

	read(fd, A, n * n * sizeof(double));
	read(fd, B, n * n * sizeof(double));
	read(fd, C, n * n * sizeof(double));
	read(fd, D, n * n * sizeof(double));

	double *ptr = L;
	for (int i = 1; i <= n; i++) {
		read(fd, ptr, i * sizeof(double));
		ptr += i;
		lseek(fd, (n - i) * sizeof(double), SEEK_CUR);
	}

	ptr = U;
	for (int i = 1; i <= n; i++) {
		read(fd, ptr, i * sizeof(double));
		ptr += i;
		lseek(fd, (n - i) * sizeof(double), SEEK_CUR);
	}

	read(fd, res, n * n * sizeof(double));

	close(fd);
}

INLINE
void common(int rank, int n, int t, double *A, double *B, double *C, double *D,
		double *L, double *U, double *aux1, double **result)
{
	const double elem = n * n;
	double sums[2];
	const int total = n * n;
	const int slice = total / t;
	double dur;
	double durrel;

	/*
	 * Distribución de datos
	 */

	MPI_Scatter(A, slice, MPI_DOUBLE, A, slice,
			MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(B, total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(C, total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(D, slice, MPI_DOUBLE, D, slice,
			MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(U, n * (n + 1) >> 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	/* El caso especial de L */
	{
		int sendcounts[t];
		int displs[t];
		for (int i = 0; i < t; i++) {
			displs[i] = L_FIL(n / t * i, 0);
			sendcounts[i] = L_FIL(n / t * (i + 1), 0) - displs[i];
		}
		MPI_Scatterv(L, sendcounts, displs, MPI_DOUBLE,
			L, sendcounts[rank], MPI_DOUBLE,
			0, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * Cómputo
	 */

	durrel = dwalltime();

	double *AB = aux1;
	multiply(AB, A, B, n, t, rank);


	double *LC = A;
	sums[0] = multiply_ll(LC, L, C, n, t, rank);

	double *DU = B;
	sums[1] = multiply_ru(DU, D, U, n, t, rank);


	double *ABpLC = A;
	add(ABpLC, AB, LC, n, t, rank);

	double *res = AB;

	dur = dwalltime() - durrel;

	MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	durrel = dwalltime();

	add(res, ABpLC, DU, n, t, rank);
	scale(res, sums[0] * sums[1] / (elem * elem), n, t, rank);

	dur += dwalltime() - durrel;
	printf("T(rank %d) = %f [s] \n", rank, dur);

	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * Gather del resultado
	 */

	MPI_Gather(res, slice, MPI_DOUBLE, res, slice,
			MPI_DOUBLE, 0, MPI_COMM_WORLD);

	*result = res;

}

INLINE
void master(const char *filename, int n, int t)
{
	double *A, *B, *C, *D, *U, *L, *R1, *result, *expected_result;
	double dur;

	A = malloc(n * n * sizeof(double));
	B = malloc(n * n * sizeof(double));
	C = malloc(n * n * sizeof(double));
	D = malloc(n * n * sizeof(double));
	L = malloc(n * (n + 1) / 2 * sizeof(double));
	U = malloc(n * (n + 1) / 2 * sizeof(double));
	R1 = malloc(n * n * sizeof(double));

	if (filename) {
		expected_result = malloc(n * n * sizeof(double));
		load(filename, A, B, C, D, L, U, expected_result, n);
	} else {
		for (int i = 0; i < n * n; i++)
			A[i] = 1;
		memcpy(B, A, n * n * sizeof(double));
		memcpy(C, A, n * n * sizeof(double));
		memcpy(D, A, n * n * sizeof(double));
		memcpy(U, A, n * (n + 1) /2 * sizeof(double));
		memcpy(L, A, n * (n + 1) /2 * sizeof(double));
	}

	dur = dwalltime();

	common(0, n, t, A, B, C, D, L, U, R1, &result);

	printf("T = %f [s] \n", dwalltime() - dur);

	if (filename)
	{
		double rms[2];
		error(expected_result, result, n, rms);
		fprintf(stderr, "RMS*RMS = %.25f\n", rms[0]);
		fprintf(stderr, "promedio = %.25f\n", rms[1]);

		free(expected_result);
	}

	free(A);
	free(B);
	free(C);
	free(D);
	free(U);
	free(L);
	free(R1);
}

INLINE
void slave(int rank, int n, int t)
{
	double *A, *B, *C, *D, *U, *L, *R1, *result;
	A = malloc(n * n / t * sizeof(double));
	B = malloc(n * n *  sizeof(double));
	C = malloc(n * n * sizeof(double));
	D = malloc(n * n / t * sizeof(double));
	L = malloc((L_FIL(n / t * (rank + 1), 0) - L_FIL(n / t * rank, 0)) * sizeof(double));
	U = malloc(n * (n + 1) / 2 * sizeof(double));
	R1 = malloc(n * n / t * sizeof(double));

	common(rank, n, t, A, B, C, D, L, U, R1, &result);

	free(A);
	free(B);
	free(C);
	free(D);
	free(U);
	free(L);
	free(R1);
}

int main(int argc, char **argv)
{
	int rank, t, n;
	const char *filename = NULL;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Comm_size(MPI_COMM_WORLD, &t);

#ifdef _OPENMP
	if (argc != 3 && argc != 4) {
		if (rank == 0)
			fprintf(stderr, "Uso:\nparalelo N hilos [archivo]\n");
		goto skip;
	}
	omp_set_num_threads(atoi(argv[2]));
	if (argc == 4)
		filename = argv[3];
#else
	if (argc != 2 && argc != 3) {
		if (rank == 0)
			fprintf(stderr, "Uso:\nparalelo N [archivo]\n");
		goto skip;
	}
	if (argc == 3)
		filename = argv[2];
#endif

	n = atoi(argv[1]);

	if (rank == 0)
		master(filename, n, t);
	else
		slave(rank, n, t);

skip:
	MPI_Finalize();
	return 0;
}
