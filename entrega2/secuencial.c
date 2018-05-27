#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

#define TAG 1
#define INLINE __attribute__((always_inline)) inline

/* Defines para acceso a matrices triangulares */
#define U_FIL(i,j) (i * n + j - (i * (i + 1)) / 2)
#define U_COL(i,j) (i + (j * (j + 1)) / 2)
#define L_FIL(i,j) (j + (i * (i + 1)) / 2)

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

/*
 * Suma A y B, dejando el resultado en C. A y B deben estar almacenadas por
 * filas. C puede ser A o B.
 */
INLINE
void add(double *C, const double *A, const double *B, int n)
{
	for (int i = 0; i < n * n; i++)
		C[i] = A[i] + B[i];
}

/*
 * Multiplica cada elemento de A por factor. Deja el resultado directamente en
 * A.
 */
INLINE
void scale(double *A, double factor, int n)
{
	for (int i = 0; i < n * n; i++)
		A[i] *= factor;
}

/*
 * Suma todos los elementos de A. Añade la suma a la variable res atómicamente,
 * usando el mútex mutex.
 */
INLINE
double sum(const double *A, int n)
{
	double s = 0;

	for (int i = 0; i < n * n; i++)
		s += A[i];

	return s;
}

INLINE
void multiply(double *C, const double *restrict A,
	      const double *restrict B, int n)
{
	double c;

	/* Multiplicación convencional fila * columna */
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			c = 0;
			for (int k = 0; k < n; k++) {
				c += A[i * n + k] * B[j * n + k];
			}
			C[i * n + j] = c;
		}
	}

}

/*
 * Calcula el producto A * B donde A es una matriz triangular inferior ordenada
 * por filas.
 * A la vez, calcula la suma de los elementos de A.
 * El resultado se almacena en C. C != B != A.
 */
INLINE
double multiply_ll(double *restrict C, const double *restrict A,
		 const double *restrict B, int n)
{
	double c;
	double s = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			c = 0;
			for (int k = 0; k <= i; k++)
				c += A[L_FIL(i, k)] * B[j * n + k];

			C[i * n + j] = c;
			if (i >= j)
				s += A[L_FIL(i,j)];
		}
	}
	return s;
}

/*
 * Calcula el producto A * B donde B es una matriz triangular superior ordenada
 * por columnas.
 * A la vez, calcula la suma de los elementos de B.
 * El resultado se almacena en C. C != B != A.
 */
INLINE
double multiply_ru(double *restrict C, const double *restrict A,
		 const double *restrict B, int n)
{
	double c;
	double s = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			c = 0;
			for (int k = 0; k <= j; k++)
				c += A[i * n + k] * B[U_COL(k, j)];

			C[i * n + j] = c;
			if (i <= j)
				s += B[U_COL(i, j)];
		}
	}

	return s;
}


/*
 * Compara las matrices A y B con una tolerancia epsilon. A y B deben estar 
 * almacenadas de la misma forma.
 */
double rms(double *A, double *B, int n)
{
	double res = 0;
	for (int i = 0; i < n * n; i++) {
		double d = (A[i] - B[i]);
		d *= d;
		res += d;
	}

	return res / (n * n);
}


int main(int argc, char **argv)
{
	bool use_file;

	/* Dimensión de matrices */
	int n;

	/* Matrices entrada */
	double *A, *B, *C, *D, *U, *L;
	/* Matriz auxiliar */
	double *R1;

	/* Puntero a resultado final */
	double *result;

	/* Matriz con el resultado dado como entrada */
	double *given_result = NULL;

	/* Tiempos */
	double ti, tf;

	if (argc != 2 && argc != 3) {
		printf("secuencial n [archivo]\n");
		return -1;
	}

	use_file = argc == 3;

	/* Dimensión de las matrices */
	n = atoi(argv[1]);

	A = malloc(sizeof(double) * n * n);
	B = malloc(sizeof(double) * n * n);
	C = malloc(sizeof(double) * n * n);
	D = malloc(sizeof(double) * n * n);
	L = malloc(sizeof(double) * (n * (n + 1) / 2));
	U = malloc(sizeof(double) * (n * (n + 1) / 2));
	R1 = malloc(sizeof(double) * n * n);

	/*
	 * Si se pasa un archivo, usar para cargar la matrices. Si no,
	 * inicializarlas en 1.
	 */
	if (use_file) {
		int fd = open(argv[2], O_RDONLY);

		given_result = malloc(sizeof(double) * n * n);

		if (fd < 1) {
			perror("load");
			exit(-1);
		}

		read(fd, A, n * n * sizeof(double));
		read(fd, B, n * n * sizeof(double));
		read(fd, C, n * n * sizeof(double));
		read(fd, D, n * n * sizeof(double));

		double *ptr = L;
		for (int i = 1; i <= n; i++) {
			read(fd, ptr, i * sizeof(double));
			ptr += i;
			lseek(fd, (n - i) * sizeof(double), SEEK_CUR);
		}

		ptr = U;
		for (int i = 1; i <= n; i++) {
			read(fd, ptr, i * sizeof(double));
			ptr += i;
			lseek(fd, (n - i) * sizeof(double), SEEK_CUR);
		}

		read(fd, given_result, n * n * sizeof(double));

		close(fd);
	} else {
		const size_t size = n * n * sizeof(double);
		const size_t tsize = (n * (n + 1)) / 2 * sizeof(double);

		for (int i = 0; i < n * n; i++)
			A[i] = 1;

		memcpy(B, A, size);
		memcpy(C, A, size);
		memcpy(D, A, size);
		memcpy(L, A, tsize);
		memcpy(U, A, tsize);
	}

	ti = dwalltime();
	{
		const double elem = n * n;
		double avg_l, avg_u;

		double *AB = R1;
		multiply(AB, A, B, n);

		double *LC = A;
		avg_l = multiply_ll(LC, L, C, n) / elem;

		double *DU = B;
		avg_u = multiply_ru(DU, D, U, n) / elem;

		double *ABpLC = A;
		add(ABpLC, AB, LC, n);

		result = ABpLC;
		add(result, ABpLC, DU, n);
		scale(result, avg_u * avg_l, n);

	}
	tf = dwalltime();

	printf("T = %f [s]\n", tf - ti);

	if (use_file)
		fprintf(stderr, "RMS*RMS = %f\n", rms(given_result, result, n));

	free(A);
	free(B);
	free(C);
	free(D);
	free(L);
	free(U);
	free(R1);
	free(given_result);

	return 0;
}
