#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <float.h>

#include <sys/time.h>


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


int main(int argc, char **argv)
{
	int n;
	double *A;
	double *B;
	double avg = 0;
	double max = -DBL_MAX;
	double min = DBL_MAX;

	if (argc != 2) {
		fprintf(stderr, "Uso:\nsecuencial N\n");
		return -1;
	}

	n = atoi(argv[1]);

	A = malloc(sizeof(*A) * n * n);
	B = malloc(sizeof(*B) * n * n);

	init(A, n);

	double dur = dwalltime();
	for (int i = 0; i < n * n; i++) {
		if (A[i] < min)
			min = A[i];
		if (A[i] > max)
			max = A[i];
		avg += A[i];
	}
	avg = avg / (n * n);

	for (int i = 0; i < n * n; i++) {
		double a = A[i];
		if (a > avg)
			B[i] = max;
		else if (a < avg)
			B[i] = min;
		else
			B[i] = avg; // Dudo que esto se ejecute
	}

	dur = dwalltime() - dur;

	fprintf(stderr, "MAX\t%f\n", max);
	fprintf(stderr, "MIN\t%f\n", min);
	fprintf(stderr, "AVG\t%f\n", avg);

	printf("T = %f [s]\n", dur);

	free(A);
	free(B);

	return 0;
}

