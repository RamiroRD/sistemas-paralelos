#!/bin/sh

make 1>&2

octave generar.m 512 1>&2
octave generar.m 1024 1>&2
octave generar.m 2048 1>&2

echo "Ejercicio 1 - MPI"
for T in 2 4 ; do
	for N in 512 1024 2048 ; do
		echo -n N = $N '\t'T = $T '\t'
		1/pthreads/producto $T $N
	done
done

echo "Ejercicio 2 - Híbrido MPI + OpenMP"
for N in 512 1024 2048 ; do
		echo -n N = $N '\t'SEQ '\t'
		1/openmp/producto-secuencial 1 $N
done

echo FIN

