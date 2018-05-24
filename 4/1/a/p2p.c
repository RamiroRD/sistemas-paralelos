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
inline void multiply(const double *A, const double *B, double *C, int n, int t)
{
	const int rowslice = n / t;
	double c;
	for (int i = 0; i < rowslice; i++) {
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
	double ti, tf;
	const int slice = n * n / t;
	const size_t size = n * n;
	MPI_Request requests[2 * (t - 1)];

	A = malloc(size * sizeof(double));
	B = malloc(size * sizeof(double));
	C = malloc(size * sizeof(double));

	one(A, n);
	one(B, n);

	ti = dwalltime();

	MPI_Request *req = requests;
	for (int dest = 1; dest < t; dest++) {
		MPI_Isend(A + dest * slice, slice, MPI_DOUBLE, dest, TAG, MPI_COMM_WORLD, req++);
		MPI_Isend(B, n * n, MPI_DOUBLE, dest, TAG, MPI_COMM_WORLD, req++);
	}

	multiply(A, B, C, n, t);

	MPI_Waitall(2 * (t - 1), requests, MPI_STATUSES_IGNORE);

	req = requests;
	for (int slave = 1; slave < t; slave++) {
		MPI_Irecv(C + slave * slice, slice, MPI_DOUBLE, slave, TAG, MPI_COMM_WORLD, req++);
	}

	MPI_Waitall(t - 1, requests, MPI_STATUSES_IGNORE);
	tf = dwalltime();

	// wait_for_gdb();

	fprintf(stderr, check(C, n) ? "OK\n" : "ERROR\n");

	printf("T = %f\n", tf - ti);
}

INLINE
inline void slave(int rank, int n, int t)
{
	double *A, *B, *C;
	const int slice = n * n / t;
	const size_t size = n * n * sizeof(double);

	A = malloc(slice * sizeof(double));
	B = malloc(size);
	C = malloc(slice * sizeof(double));

	MPI_Recv(A, slice, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(B, n * n, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	multiply(A, B, C, n, t);

	MPI_Send(C, slice, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
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
