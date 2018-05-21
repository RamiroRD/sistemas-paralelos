#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <float.h>

#include <unistd.h>
#include <sys/time.h>

#include <mpi.h>

#define INLINE __attribute__((always_inline))

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;

	return sec;
}

void init(double *V, uint64_t n)
{
	for (uint64_t i = 0; i < n * n; i++)
		V[i] = 1000.0 - 1000.0 / (n * n) * i;
}

INLINE
inline void common(int rank, int n, int t, double *A, double *B)
{
	double max, min, sum, avg;
	double lmax = DBL_MIN;
	double lmin = DBL_MAX;
	double lsum = 0;

	for (int i = 0; i < n * n / t; i++) {
		if (A[i] < lmin)
			lmin = A[i];
		if (A[i] > lmax)
			lmax = A[i];
		lsum += A[i];
	}

	MPI_Allreduce(A, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(A, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(A, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	avg = sum / (n * n);

	for (int i = 0; i < n * n / t; i++) {
		double a = A[i];
		if (a > avg)
			B[i] = max;
		else if (a < avg)
			B[i] = min;
		else
			B[i] = avg; // Dudo que esto se ejecute
	}

}


INLINE inline void master(int n, int t)
{
	double *A, *B;
	double dur;

	A = malloc(sizeof(double) * n * n);
	B = malloc(sizeof(double) * n * n);

	init(A, n);

	dur = dwalltime();
	MPI_Scatter(A, n * n / t, MPI_DOUBLE, A, n * n / t, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	common(0, n, t, A, B);

	dur = dwalltime() - dur;

	printf("T = %f [s]\n", dur);

	free(A);
	free(B);
}

INLINE inline void slave(int rank, int n, int t)
{
	double *A, *B;

	A = malloc(sizeof(double) * n * n / t);
	B = malloc(sizeof(double) * n * n / t);

	MPI_Scatter(NULL, 0, MPI_DOUBLE, A, n * n / t, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	common(rank, n, t, A, B);

	free(A);
	free(B);
}

int main(int argc, char **argv)
{
	int rank, t, n;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Comm_size(MPI_COMM_WORLD, &t);

	if (argc != 2) {
		if (rank == 0)
			fprintf(stderr, "Uso:\nparalelo <N>\n");
		goto skip;
	}

	n = atoi(argv[1]);

	if (rank == 0)
		master(n, t);
	else
		slave(rank, n, t);

      skip:
	MPI_Finalize();
	return 0;
}
