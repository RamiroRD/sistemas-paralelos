#include <stdio.h>
#include <mpi.h>

#define TAG 1
#define INLINE __attribute__((always_inline))

INLINE
inline void master(int t)
{
	for(int i = 1; i < t; i++) {
		int num = 2 * i;
		printf("MAESTRO enviando a %d.\n", i);
		MPI_Send(&num, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
	}
}

INLINE
inline void slave(int rank, int t) 
{
		int num;
		MPI_Recv(&num, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("ESCLAVO %i recibe %i.\n", rank, num);
}


int main(int argc, char **argv)
{
	int rank, t;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Comm_size(MPI_COMM_WORLD, &t);

	if (rank == 0)
		master(t);
	else
		slave(rank, t);

	MPI_Finalize();

	return 0;
}
