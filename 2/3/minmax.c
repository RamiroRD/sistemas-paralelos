#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

/* Globales */
double *V;
int n, t;

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
	printf("minmax <n> <t>\n");
}

bool ispwtwo(int n)
{
	while ((n & 1) == 0) {
		n >>= 1;
	}
	return n == 1;
}

struct slave_args {
	int base;
	double min, max;
};

void *slave_routine(void *args)
{
	struct slave_args *arg = (struct slave_args*) args;
	double min = INFINITY, max = -INFINITY;
	int base = arg->base;

	for (int i = base; i < base + n / t; i++) {
		if (V[i] < min)
			min = V[i];
		if (V[i] > max)
			max = V[i];
	}

	arg->min = min;
	arg->max = max;
	return args;
}

int main(int argc, char **argv)
{	
	pthread_t *threads;
	struct slave_args *args;

	if (argc < 3) {
		usage();
		return -1;
	}

	n = atoi(argv[1]);
	t = atoi(argv[2]);

	if (!ispwtwo(t) || !ispwtwo(n)) {
		printf("Sólo potencias de dos.\n");
		return -1;
	}

	V = (double*) malloc(sizeof(double) * n);

	srand(time(NULL));
	for (int i = 0; i < n; i++)
		V[i] = rand();

	threads = (pthread_t*) malloc(sizeof(pthread_t) * t);
	args = (struct slave_args*) malloc(sizeof(struct slave_args) * t);
	
	double dur = dwalltime();
	for (int i = 0; i < t; i++) {
		args[i].base = i * n / t;
		pthread_create(threads + i, NULL, slave_routine, args + i);
	}

	for (int i = 0; i < t; i++)
		pthread_join(threads[i], NULL);
	dur = dwalltime() - dur;

	printf("%f\n", dur);
	
	free(threads);
	free(args);
	free(V);

	return 0;
}

