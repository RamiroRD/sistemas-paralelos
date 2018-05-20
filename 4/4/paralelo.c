#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

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

void merge(double *V, double *buffer, uint64_t al, const uint64_t ar, uint64_t bl,
	 const uint64_t br)
{
	const uint64_t base = al;
	uint64_t res = 0;
	while (al <= ar && bl <= br)
		buffer[res++] = (V[al] < V[bl]) ? V[al++] : V[bl++];
	while (al <= ar)
		buffer[res++] = V[al++];
	while (bl <= br)
		buffer[res++] = V[bl++];
	memcpy(V + base, buffer, sizeof(double) * (br - base + 1));
}

void do_sort(double *V, double *buffer, uint64_t l, uint64_t r)
{
	if (r - l == 0)
		return;
	do_sort(V, buffer, l, l + (r - l) / 2);
	do_sort(V, buffer, l + (r - l) / 2 + 1, r);
	merge(V, buffer, l, l + (r - l) / 2, l + (r - l) / 2 + 1, r);
}

void init(double *V, uint64_t n)
{
	for (uint64_t i = 0; i < n; i++)
		V[i] = 1000.0 - 1000.0 / n * i;
}

bool check(double *V, uint64_t n)
{
	for (uint64_t i = 0; i < n - 1; i++)
		if (V[i] > V[i + 1])
			return false;
	return true;
}

void wait_for_gdb()
{
	int edit_me = 0;
	while (!edit_me)
		sleep(1);
}

INLINE
inline void common(int n, int t, int rank, double *V, double *merge_buffer)
{
	int step = 1;
	while(rank % step == 0 && step <= t) {
		if (step == 1) {
			do_sort(V, merge_buffer, 0, n / t - 1);
		} else {
			/* Recibimos poniendo todo al final de V */
			uint64_t local_slice = n / t * (step >> 1);
			/* Usamos el tag para diferenciar los sends */
			MPI_Recv(V + local_slice, local_slice, MPI_DOUBLE, MPI_ANY_SOURCE, step >> 1, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
			
			merge(V, merge_buffer, 0, local_slice - 1, local_slice, 2 * local_slice - 1);
		}

		if ((rank % step == 0) && ((rank / step) % 2 != 0))
			MPI_Send(V, n / t * step, MPI_DOUBLE, rank - step, step, MPI_COMM_WORLD);
		step <<= 1;
	}
}


INLINE
inline void master(int n, int t)
{
	double *V;
	double *merge_buffer;
	double dur;

	V = malloc(sizeof(*V) * n);
	merge_buffer = malloc(sizeof(*V) * n);

	init(V, n);

	dur = dwalltime();
	MPI_Scatter(V, n / t, MPI_DOUBLE, V, n / t, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	common(n, t, 0, V, merge_buffer);

	fprintf(stderr, check(V, n) ? "OK\n" : "ERROR\n");
	printf("T = %f [s]\n", dwalltime() - dur);
}

/* TODO: Optimizar memoria */
/* TODO: Microoptimizar */
INLINE
inline void slave(int rank, int n, int t) 
{
	double *V;
	double *merge_buffer;

	V = malloc(sizeof(*V) * n);
	merge_buffer = malloc(sizeof(*V) * n);

	MPI_Scatter(NULL, 0, MPI_DOUBLE, V, n / t , MPI_DOUBLE, 0, MPI_COMM_WORLD);

	common(n, t, rank, V, merge_buffer);

	free(V);
	free(merge_buffer);
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

