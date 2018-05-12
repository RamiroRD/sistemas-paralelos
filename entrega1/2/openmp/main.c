#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Defines para acceso a matrices triangulares */
#define U_FIL(i,j) (i * n + j - (i * (i + 1)) / 2)
#define U_COL(i,j) (i + (j * (j + 1)) / 2)
#define L_FIL(i,j) (j + (i * (i + 1)) / 2)

bool use_file;

/* Cantidad de hilos y dimensión de matrices */
int t, n;

/* Matrices entrada y alguna de sus transpuestas */
double *A, *B, *C, *D, *E, *F, *U, *L;
double *AT, *BT, *CT, *ET, *FT, *UT;

/* Matrices intermedias. No se reservan memoria exclusiva para estos*/
double *AA, *AAC;
double *LB, *LBE;
double *DU, *DUF;
double *LBEpDUF;

/* Puntero a resultado final */
double *result;

/* Matriz con el resultado dado como entrada */
double *given_result;

/* Sumas de U, L y B */
double sum_u = 0, sum_l = 0, sum_b = 0;

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

/*
 * Traspone src, dejando el resultado en dst. src y dst deben ser distintos.
 */
void transpose(double *restrict dst, const double *restrict src, int n)
{
#pragma omp parallel for
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			dst[j * n + i] = src[i * n + j];
}


void transpose_upper(double *restrict dst, const double *restrict src,
		     int n)
{
#pragma omp parallel for schedule(dynamic,4)
	for (int j = 0; j < n; j++)
		for (int i = 0; i <= j; i++)
			dst[U_COL(i, j)] = src[U_FIL(i, j)];

}

/*
 * Suma A y B, dejando el resultado en C. A y B deben estar almacenadas por
 * filas. C puede ser A o B.
 */
void add(double *C, const double *A, const double *B, int n)
{
#pragma omp parallel for
	for (int i = 0; i < n * n; i++)
		C[i] = A[i] + B[i];
}

/*
 * Multiplica cada elemento de A por factor. Deja el resultado directamente en
 * A.
 */
void scale(double *A, double factor, int dim)
{
#pragma omp parallel for
	for (int i = 0; i < dim; i++)
		A[i] *= factor;
}

/*
 * Suma todos los elementos de A.
 */
void sum(const double *A, double *res, int dim)
{
	double partial = 0;
#pragma omp parallel for reduction(+:partial)
	for (int i = 0; i < dim; i++)
		partial += A[i];

	*res = partial;
}

void multiply(double *C, const double *restrict A,
	      const double *restrict B, int n)
{
	/* Multiplicación convencional fila * columna */
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			double partial = 0;
			for (int k = 0; k < n; k++)
				partial += A[i * n + k] * B[j * n + k];
			C[i * n + j] = partial;
		}
}

/*
 * Calcula el producto A * B donde A es una matriz triangular inferior ordenada
 * por filas.
 * El resultado se almacena en C. C != B != A.
 */
void multiply_ll(double *restrict C, const double *restrict A,
		 const double *restrict B, int n)
{
	/* Multiplicación convencional fila * columna */
#pragma omp parallel for schedule(dynamic,4)
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double partial = 0;
			/* Cuando k > i los elementos de A valen 0, por lo tanto se deja de sumar. */
			for (int k = 0; k <= i; k++)
				partial += A[L_FIL(i, k)] * B[j * n + k];
			C[i * n + j] = partial;
		}
	}
}

/*
 * Calcula el producto A * B donde B es una matriz triangular superior ordenada
 * por columnas.
 * El resultado se almacena en C. C != B != A.
 */
void multiply_ru(double *restrict C, const double *restrict A,
		 const double *restrict B, int n)
{
	/* Multiplicación convencional fila * columna */
#pragma omp parallel for schedule(dynamic,4)
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			double partial = 0;
			/* Cuando k > j los elementos de B valen 0, por lo tanto se deja de sumar. */
			for (int k = 0; k <= j; k++)
				partial += A[i * n + k] * B[U_COL(k, j)];
			C[i * n + j] = partial;
		}
}


void load(const char *file)
{
	int fd = open(file, O_RDONLY);
	if (fd < 1) {
		perror("load");
		exit(-1);
	}
	read(fd, A, n * n * sizeof(double));
	read(fd, B, n * n * sizeof(double));
	read(fd, C, n * n * sizeof(double));
	read(fd, D, n * n * sizeof(double));
	read(fd, E, n * n * sizeof(double));
	read(fd, F, n * n * sizeof(double));

	double *ptr = L;
	for (int i = 1; i <= n; i++) {
		read(fd, ptr, i * sizeof(double));
		ptr += i;
		lseek(fd, (n - i) * sizeof(double), SEEK_CUR);
	}

	ptr = U;
	for (int i = n; i >= 1; i--) {
		lseek(fd, (n - i) * sizeof(double), SEEK_CUR);
		read(fd, ptr, i * sizeof(double));
		ptr += i;
	}

	read(fd, given_result, n * n * sizeof(double));

	close(fd);
}

/*
 * Inicializa todas las matrices en 1.
 */
void ones()
{
	const size_t size = n * n * sizeof(double);
	const size_t tsize = (n * (n + 1)) / 2 * sizeof(double);

	for (int i = 0; i < n * n; i++)
		A[i] = 1;

	memcpy(B, A, size);
	memcpy(C, A, size);
	memcpy(D, A, size);
	memcpy(E, A, size);
	memcpy(F, A, size);
	memcpy(L, A, tsize);
	memcpy(U, A, tsize);
}

/*
 * Compara las matrices A y B con una tolerancia epsilon. A y B deben estar
 * almacenadas de la misma forma.
 */
bool compare(double *A, double *B, int n, double epsilon)
{
	for (int i = 0; i < n * n; i++) {
		if (fabs(A[i] - B[i]) > (A[i] * epsilon))
			return false;
	}

	return true;
}

int main(int argc, char **argv)
{
	/* Tiempos */
	double ti, tf;

	if (argc != 3 && argc != 4) {
		printf("producto T n [archivo]\n");
		return -1;
	}

	use_file = argc == 4;
	/* Cantidad de hilos */
	t = atoi(argv[1]);
	/* Dimensión de bloque (en elementos) */
	n = atoi(argv[2]);

	/*
	 * Reservamos memoria para las matrices A, B, C, D, E, F, L, U y para
	 * la representación ordenada por columnas de A, B, C, E, F, U.
	 */
	A = malloc(sizeof(double) * n * n);
	B = malloc(sizeof(double) * n * n);
	C = malloc(sizeof(double) * n * n);
	D = malloc(sizeof(double) * n * n);
	E = malloc(sizeof(double) * n * n);
	F = malloc(sizeof(double) * n * n);
	L = malloc(sizeof(double) * (n * (n + 1) / 2));
	U = malloc(sizeof(double) * (n * (n + 1) / 2));
	AT = malloc(sizeof(double) * n * n);
	BT = malloc(sizeof(double) * n * n);
	CT = malloc(sizeof(double) * n * n);
	ET = malloc(sizeof(double) * n * n);
	FT = malloc(sizeof(double) * n * n);
	UT = malloc(sizeof(double) * (n * (n + 1) / 2));
	given_result = malloc(sizeof(double) * n * n);

	/*
	 * Si se pasa un archivo, usar para cargar la matrices. Si no,
	 * inicializarlas en 1.
	 */
	if (use_file)
		load(argv[3]);
	else
		ones();


	/*
	 * Seteamos la cantidad de hilos.
	 */
#ifdef _OPENMP
	omp_set_num_threads(t);
#endif

	ti = dwalltime();
	/* Promedio de u */
	sum(U, &sum_u, (n * (n + 1)) / 2);
	/* Promedio de l */
	sum(L, &sum_l, (n * (n + 1)) / 2);
	/* Promedio de b */
	sum(B, &sum_b, n * n);

	double avg_u = sum_u / (n * n);
	double avg_l = sum_l / (n * n);
	double avg_b = sum_b / (n * n);

	/* Transpuestas */
	transpose(AT, A, n);
	transpose(BT, B, n);
	transpose(CT, C, n);
	transpose(ET, E, n);
	transpose(FT, F, n);
	transpose_upper(UT, U, n);

	/*  AA */
	AA = C;			/* Reutilizamos el espacio de C */
	multiply(AA, A, AT, n);

	/* AAC */
	AAC = A;		/* Reutilizamos A */
	multiply(AAC, AA, CT, n);

	/* ulAAC */
	scale(AAC, avg_u * avg_l, n * n);

	/* LB */
	LB = C;			/* Reutilizamos el espacio de C de vuelta */
	multiply_ll(LB, L, BT, n);

	/* LBE */
	LBE = B;		/* Reutilizamos el espacio de B */
	multiply(LBE, LB, ET, n);

	/* DU */
	DU = C;			/* Reutilizamos el espacio de C */
	multiply_ru(DU, D, UT, n);

	/* DUF */
	DUF = E;		/* Reutilizamos el espacio de E */
	multiply(DUF, DU, FT, n);

	/* b/(LBE + DUF) en el espacio de C */
	LBEpDUF = C;
	add(LBEpDUF, LBE, DUF, n);
	scale(LBEpDUF, avg_b, n);

	result = A;
	/* Resultado final en A */
	add(result, AAC, LBEpDUF, n);

	tf = dwalltime();

	printf("T = %f [s]\n", tf - ti);

	if (use_file) {
		if (compare(given_result, result, n, 0.01))
			fprintf(stderr, "OK\n");
		else
			fprintf(stderr, "ERROR\n");
	}

	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	free(L);
	free(U);
	free(AT);
	free(BT);
	free(CT);
	free(ET);
	free(FT);
	free(UT);
	free(given_result);

	return 0;
}
