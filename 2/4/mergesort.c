#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>

/* Globales */
int n, t;
int *V;

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
	printf("mergesort <n> <t>\n");
}

bool ispwtwo(int n)
{
	while ((n & 1) == 0) {
		n >>= 1;
	}
	return n == 1;
}


void merge(int *V, int *buffer, int al, const int ar, int bl,
	 const int br)
{
	const int base = al;
	int res = 0;
	while (al <= ar && bl <= br)
		buffer[res++] = (V[al] < V[bl]) ? V[al++] : V[bl++];
	while (al <= ar)
		buffer[res++] = V[al++];
	while (bl <= br)
		buffer[res++] = V[bl++];
	memcpy(V + base, buffer, sizeof(int) * (br - base + 1));
}

void do_sort(int *V, int *buffer, int l, int r)
{
	if (r - l == 0)
		return;
	do_sort(V, buffer, l, l + (r - l) / 2);
	do_sort(V, buffer, l + (r - l) / 2 + 1, r);
	merge(V, buffer, l, l + (r - l) / 2, l + (r - l) / 2 + 1, r);
}

void print(int *V, int n) {
	for (int i = 0; i < n; i++)
		printf("%d\t", V[i]);

	printf("\r\n");
}

void init(int *V, int n)
{
	srand(time(NULL));
	for (int i = 0; i < n; i++)
		V[i] = rand() % 100;
}

bool check(int *V, int n)
{
	for (int i = 0; i < n - 1; i++)
		if (V[i] > V[i + 1])
			return false;
	return true;
}

struct worker_args {
	int id;
	int t;
};

void *worker(void *args)
{
	struct worker_args* arg = (struct worker_args*) args;
	const int id = arg->id;
	const int slice = n / arg->t;
	int *buffer = (int*) malloc(sizeof(int) * slice);
	do_sort(V, buffer, id * slice, (id + 1) * slice - 1);


	free(buffer);
	return args;
}

void *merger(void *args)
{
	struct worker_args* arg = (struct worker_args*) args;
	const int id = arg->id;
	const int slice = n / arg->t;
	int *buffer = (int*) malloc(sizeof(int) * slice);
	merge(V, buffer, id * slice, id * slice + slice / 2 - 1, id * slice + slice / 2,  (id + 1) * slice - 1);
	free(buffer);
	return NULL;
}

int main(int argc, char **argv)
{
	pthread_t *threads;
	struct worker_args *args;

	if (argc < 3) {
		usage();
		return -1;
	}

	n = atoi(argv[1]);
	t = atoi(argv[2]);

	if (!ispwtwo(t) || !ispwtwo(n)) {
		printf("Solo potencias de dos.\n");
		return -1;
	}

	V = (int*) malloc(sizeof(int) * n);
	threads = (pthread_t*) malloc(sizeof(pthread_t) * t);
	args = (struct worker_args*) malloc(sizeof(struct worker_args) * t);

	init(V, n);

	double dur = dwalltime();
	for (int i = 0; i < t; i++) {
		args[i].id = i;
		args[i].t = t;
		pthread_create(threads + i, NULL, worker, args + i);
	}

	for (int i = 0; i < t; i++)
		pthread_join(threads[i], NULL);

	for (int i = t / 2; i > 1; i >>= 1) {
		for (int j = 0; j < i; j++) {
			args[j].t = i;
			pthread_create(threads + j, NULL, merger, args + j);
		}
		for (int j = 0; j < i; j++)
			pthread_join(threads[j], NULL);
	}
	
	/* Es necesario un último merge? */
	if (t != 1) {
		args[0].t = 1;
		merger(args + 0);
	}

	dur = dwalltime() - dur;

	bool ok = check(V, n);
	fprintf(stderr, ok ? "OK\n" : "ERROR\n");
	if (ok)
		printf("%f\n", dur);

	free(threads);
	free(args);
	free(V);
	return 0;
}

