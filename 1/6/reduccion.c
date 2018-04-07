#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

int main(int argc, char **argv)
{
	int p, m, n;
	double start, end;
	double *V;

	if (argc < 2) {
		fprintf(stderr, "Uso:\nreduccion <potencia>\n");
		exit(-1);
	}
	p = atoi(argv[1]);

	if (p < 1) {
		fprintf(stderr, "Potencia debe ser mayor o igual a 1.\n");
		exit(-1);
	}

	n = 1 << p;
	V = (double*) (malloc(sizeof(double) * n));	
	for (int i = 0; i < n; i++)
		V[i] = (double) i + 1;	

	start = dwalltime();
	m = 1;	
	while (n > 1) {
		for(int i = 0; i < n; i++)
			V[2 * i * m] /= V[(2 * i + 1) * m];
		m <<= 1;
		n >>= 1;
	}
	end = dwalltime();

	printf("Resultado: %f\n", V[0]);
	printf("Tiempo en segundos: %f\n", end - start);

	free(V);

	return 0;
}

