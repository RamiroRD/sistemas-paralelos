#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

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

/* Puntero a resultado final */
double *result;

/* Matriz con el resultado dado como entrada */
double *given_result;

/* Matrices intermedias. No se reservan memoria exclusiva para estos. */
double *AA, *AAC;
double *LB, *LBE;
double *DU, *DUF;

/* Sumas de U, L y B junto a sus mutexes*/
double up, lp, bp;
pthread_mutex_t up_mutex, lp_mutex, bp_mutex;
pthread_barrier_t barrier;

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
void transpose(double * restrict dst, const double * restrict src, int n, int t, int id)
{
	int slice = n / t;
	for (int i = id * slice; i < slice * (id + 1); i++)
		for (int j = 0; j < n; j++)
			dst[i * n + j] = src[j * n + i];
}

void multiply(double *C, const double * restrict B, const double * restrict A, int n, int t, int id)
{
	int slice = n / t;

	/* Inicializamos en cero */
	memset(C + id * slice, 0, sizeof(double) * slice * n);

	/* Multiplicación convencional fila * columna */
	for (int i = id * slice; i < slice * (id + 1); i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[j * n + k];

}

void transpose_upper(double * restrict dst, const double * restrict src, int n, int t, int id)
{
	/* Cantidad de elemtos que corresponden al hilo: slice = total_elem / t */
	int slice = n / t;
	for (int i = id * slice; i < slice * (id + 1); i++)
		for (int j = i; j < n; j++)
			dst[U_COL(i, j)] =  src[U_FIL(i, j)];
			
}

/*
 * Suma A y B, dejando el resultado en C. A y B deben estar almacenadas por
 * filas. C puede ser A o B.
 */
void add(double *C, const double *A, const double *B, int n, int t, int id)
{
	for (int i = id * (n * n) / t; i < (id + 1) * n * n / t; i++)
		C[i] = A[i] + B[i];
}

/*
 * Multiplica cada elemento de A por factor. Deja el resultado directamente en
 * A.
 */
void scale(double *A, double factor, int dim, int t, int id)
{
	for (int i = id * dim / t; i < (id + 1) * dim / t; i++)
		A[i] *= factor;
}

/*
 * Suma todos los elementos de A. Añade la suma a la variable res atómicamente,
 * usando el mútex mutex.
 */
void sum(const double *A, double *res, pthread_mutex_t *mutex, int n, int id)
{
	double partial = 0;

	for (int i = id * (n * n); i < (id + 1) * n * n; i++)
		partial += A[i];
	pthread_mutex_lock(mutex);
	*res += partial;
	pthread_mutex_unlock(mutex);
}

/*
 * Calcula el producto A * B donde A es una matriz triangular inferior ordenada
 * por filas.
 * El resultado se almacena en C. C != B != A.
 */
void multiply_ll(double * restrict C, const double * restrict A, const double * restrict B, int n, int t, int id)
{
	int slice = n / t;

	/* Multiplicación convencional fila * columna */
	for (int i = id * slice; i < slice * (id + 1); i++) {
		for (int j = 0; j < n; j++) {
			/* Cuando k > i los elementos de A valen 0, por lo tanto se deja de sumar. */
			double c = 0;
			for (int k = 0; k <= i; k++) {
				c += A[L_FIL(i,k)] * B[j * n + k];
			}
			C[i * n + j] += c;
		}
	}
}

/*
 * Calcula el producto A * B donde B es una matriz triangular superior ordenada
 * por columnas.
 * El resultado se almacena en C. C != B != A.
 */
void multiply_ru(double * restrict C, const double * restrict A, const double * restrict B, int n, int t, int id)
{
	int slice = n / t;

	/* Multiplicación convencional fila * columna */
	for (int i = id * slice; i < slice * (id + 1); i++) {
		for (int j = 0; j < n; j++) {
			double c = 0;
			/* Cuando k > j los elementos de B valen 0, por lo tanto se deja de sumar. */
			for (int k = 0; k <= j; k++) {
				c += A[i * n + k] * B[U_COL(k,j)];
			}
			C[i * n + j] = c;
		}
	}
}

void *worker(void *idp)
{
	int id = *(int*) idp;

	/* Promedio de u */
	sum(U, &up, &up_mutex, n, id);
	pthread_barrier_wait(&barrier);
	/* Promedio de l */
	sum(L, &lp, &lp_mutex, n, id);
	pthread_barrier_wait(&barrier);
	/* Promedio de b */
	sum(L, &bp, &lp_mutex, n, id);
	pthread_barrier_wait(&barrier);

	/* Transpuestas */
	transpose(AT, A, n, t, id);
	transpose(BT, B, n, t, id);
	transpose(CT, C, n, t, id);
	transpose(ET, E, n, t, id);
	transpose(FT, F, n, t, id);
	transpose_upper(UT, U, n, t, id);
	pthread_barrier_wait(&barrier);

	/*  AA */
	AA = C; /* Reutilizamos el espacio de C */
	multiply(C, A, AT, n, t, id);
	pthread_barrier_wait(&barrier);

	/* AAC */
	AAC = A; /* Reutilizamos A */
	multiply(AAC, AA, CT, n, t, id);

	/* ulAAC */
	scale(AAC, (up + lp) / (n * n * n * n), n * n, t, id);
	pthread_barrier_wait(&barrier);

	/* LB */
	LB = C;  /* Reutilizamos el espacio de C de vuelta */
	multiply_ll(LB, L, BT, n, t, id);
	pthread_barrier_wait(&barrier);

	/* LBE */
	LBE = B; /* Reutilizamos el espacio de B */
	multiply(LBE, LB, ET, n, t, id);
	pthread_barrier_wait(&barrier);

	/* DU */
	DU = C; /* Reutilizamos el espacio de C */
	multiply_ru(DU, D, UT, n, t, id);
	pthread_barrier_wait(&barrier);

	/* DUF */
	DUF = E; /* Reutilizamos el espacio de C */
	multiply(DUF, DU, FT, n, t, id);
	pthread_barrier_wait(&barrier);

	/* b/(LBE + DUF) en el espacio de C */
	add(C, LBE, DUF, n, t, id);
	scale(C, bp / (n * n) , n, t, id);

	/* Resultado final en A */
	add(A, AAC, C, n, t, id);
	result = A;


	return NULL;
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
 * Compara las matrices A y B con una tolerancia epsilon. A y B deben estar 
 * almacenadas de la misma forma.
 */
bool compare(double *A, double *B, int n, double epsilon)
{
	for (int i = 0; i < n * n; i++) {
		if (fabs(A - B) > epsilon)
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


	int ids[t];
	pthread_t threads[t];

	for (int i = 0; i < t; i++)
		ids[i] = i;

	/*
	 * Inicializamos los mútex.
	 */
	pthread_mutex_init(&up_mutex, NULL);
	pthread_mutex_init(&lp_mutex, NULL);
	pthread_mutex_init(&bp_mutex, NULL);

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

	if (use_file)
		load(argv[3]);
	/* Si no hay archivo, operamos con basura */

	/*
	 * Inicializamos la barrera y arrancamos los hilos.
	 */
	pthread_barrier_init(&barrier, NULL, t);
	ti = dwalltime();
	for (int i = 0; i < t; i++)
		pthread_create(threads + i, NULL, worker, ids + i);
	for (int i = 0; i < t; i++)
		pthread_join(threads[i], NULL);
	tf = dwalltime();

	printf("T = %f [s]\n", tf - ti);

	if (use_file) {
		if (compare(given_result, result, n, 0.1))
			printf("OK\n");
		else
			printf("ERROR\n");
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
	pthread_barrier_destroy(&barrier);
	pthread_mutex_destroy(&up_mutex);
	pthread_mutex_destroy(&lp_mutex);
	pthread_mutex_destroy(&bp_mutex);

	return 0;
}

