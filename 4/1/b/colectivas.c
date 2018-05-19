#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <unistd.h>
#include <sys/time.h>

#include <mpi.h>

#define TAG 1
#define INLINE __attribute__((always_inline))

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;

	return sec;
}

void one(double *A, int n)
{
	for (int i = 0; i < n * n; i++)
		A[i] = 1;
}

bool check(double *A, int n)
{
	bool ok = true;
	for (int i = 0; i < n * n; i++)
		ok = ok && (A[i] == n);
		
	return ok;
}

INLINE
inline void multiply(const double *A, const double *B, double *C, int n, int t, int rank)
{
	const int slice = n / t;
	double c;
	for (int i = slice * rank; i < slice * (rank + 1); i++) {
		for (int j = 0; j < n; j++) {
			c = 0;
			for (int k = 0; k < n; k++) {
				c += A[i * n + j] * B[j * n + k];
			}
			C[i * n + j] = c;
		}
	}
}

void wait_for_gdb()
{
	int edit_me = 0;
	while (!edit_me)
		sleep(1);
}


INLINE
inline void master(int n, int t)
{
	double *A, *B, *C;
	const size_t size = n * n * sizeof(double);
	double ti, tf;

	A = malloc(size);
	B = malloc(size);
	C = malloc(size);

	one(A, n);
	one(B, n);

	ti = dwalltime();

	MPI_Bcast(A, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(B, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	multiply(A, B, C, n, t, 0);

	MPI_Gather(C + 0, n * n / t, MPI_DOUBLE, C, n * n / t, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	tf = dwalltime();

	fprintf(stderr, check(C, n) ? "OK\n" : "ERROR\n");

	printf("T = %f\n", tf - ti);
	// wait_for_gdb();

}

INLINE
inline void slave(int rank, int n, int t) 
{
	double *A, *B, *C;
	const size_t size = n * n * sizeof(double);

	A = malloc(size);
	B = malloc(size);
	C = malloc(size);

	MPI_Bcast(A, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(B, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	multiply(A, B, C, n, t, rank);

	MPI_Gather(C + n * n / t * rank, n * n / t, MPI_DOUBLE, C, n * n / t, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
  
int main(int argc, char **argv)
{
	int rank, t, n;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Comm_size(MPI_COMM_WORLD, &t);

	if (argc != 2) {
		if (rank == 0)
			fprintf(stderr, "Uso:\np2p <N>\n");
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

