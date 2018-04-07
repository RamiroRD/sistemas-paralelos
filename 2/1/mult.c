#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/time.h>

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}


void usage()
{
	printf("mult <n> <t>\n");
}

bool ispwtwo(int n)
{
	while ((n & 1) == 0) {
		n >>= 1;
	}
	return n == 1;
}

struct slave_args {
	double *A;
	double *B;
	double *C;
	int n;
	int first_row;
	int rows;
};

void *slave_routine(void *args)
{
	struct slave_args *arg = (struct slave_args*) args;
	double *A = arg->A;
	double *B = arg->B;
	double *C = arg->C;
	int n = arg->n;
	for (int i = arg->first_row; i < arg->first_row + arg->rows; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[j * n + k];
	return NULL;
}

int main(int argc, char **argv)
{	
	int n, t;
	pthread_t *threads;
	struct slave_args *args;
	double *A, *B, *C;

	if (argc < 3) {
		usage();
		return -1;
	}

	n = atoi(argv[1]);
	t = atoi(argv[2]);

	if (!ispwtwo(t) || !ispwtwo(n)) {
		printf("Threads y la dimensión de la matriz deben potencia de dos.\n");
		return -1;
	}

	A = (double*) malloc(sizeof(double)*n*n);
	B = (double*) malloc(sizeof(double)*n*n);
	C = (double*) malloc(sizeof(double)*n*n);

	for (int i = 0; i < n * n; i++) {
		A[i] = B[i] = 1.0;
		C[i] = 0.0;
	}

	threads = (pthread_t*) malloc(sizeof(pthread_t) * t);
	args = (struct slave_args*) malloc(sizeof(struct slave_args) * t);
	
	double dur = dwalltime();
	for (int i = 0; i < t; i++) {
		args[i].A = A;
		args[i].B = B;
		args[i].C = C;
		args[i].n = n;
		args[i].first_row = i * n / t; 
		args[i].rows = n / t; 
		pthread_create(threads + i, NULL, slave_routine, args + i);
	}

	for (int i = 0; i < t; i++)
		pthread_join(threads[i], NULL);
	dur = dwalltime() - dur;

	int ok = 1;
	for (int i = 0; i < n * n; i++)
		ok = ok && (C[i] == n);

	fprintf(stderr, ok ? "OK\n" : "ERROR\n");
	if (ok)
		printf("%f\n", dur);
	
	free(threads);
	free(args);
	free(A);
	free(B);
	free(C);

	return 0;
}

