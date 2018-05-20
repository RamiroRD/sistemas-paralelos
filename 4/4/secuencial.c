#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include <sys/time.h>


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

int main(int argc, char **argv)
{
	uint64_t n;
	double *V;
	double *buffer;

	if (argc != 2) {
		fprintf(stderr, "Uso:\nsecuencial N\n");
		return -1;
	}

	n = atoll(argv[1]);

	V = malloc(sizeof(*V) * n);
	buffer = malloc(sizeof(*buffer) * n);

	init(V, n);

	double dur = dwalltime();
	do_sort(V, buffer, 0, n - 1);
	dur = dwalltime() - dur;

	bool ok = check(V, n);
	fprintf(stderr, ok ? "OK\n" : "ERROR\n");
	printf("T = %f [s]\n", dur);

	free(V);
	free(buffer);
	return 0;
}

