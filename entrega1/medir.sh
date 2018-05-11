#!/bin/sh

cd 1/pthreads
make 1>&2
cd ../openmp
make 1>&2

cd ../../2/
octave generar.m 512 1>&2
octave generar.m 1024 1>&2
octave generar.m 2048 1>&2
cd pthreads
make 1>&2
cd ../openmp
make 1>&2

cd ../../3/
make 1>&2
cd ../

echo "Ejercicio 1 - POSIX Threads"
for T in 2 4 ; do
	for N in 512 1024 2048 ; do
		echo -n N = $N '\t'T = $T '\t'
		1/pthreads/producto $T $N
	done
done

echo "Ejercicio 1 - OpenMP"
for N in 512 1024 2048 ; do
		echo -n N = $N '\t'SEQ '\t'
		1/openmp/producto-secuencial 1 $N
done

for T in 2 4 ; do
	for N in 512 1024 2048 ; do
		echo -n N = $N '\t'T = $T '\t'
		1/openmp/producto $T $N
	done
done

echo "Ejercicio 2 - POSIX Threads"
for T in 2 4 ; do
	for N in 512 1024 2048 ; do
		echo -n N = $N '\t'T = $T '\t'
		2/pthreads/main $T $N 2/data-$N.bin
	done
done

echo "Ejercicio 2 - OpenMP"
for N in 512 1024 2048 ; do
	echo -n N = $N '\t'SEQ '\t'
	2/openmp/main-secuencial 1 $N 2/data-$N.bin
done

for T in 2 4 ; do
	for N in 512 1024 2048 ; do
		echo -n N = $N '\t'T = $T '\t'
		2/openmp/main $T $N 2/data-$N.bin
		echo -n
	done
done


echo "Ejercicio 3 - OpenMP"
for N in 200000 2000000 20000000 200000000 2000000000 ; do
	echo -n N = $N '\t'SEQ '\t'
	3/pares-secuencial 1 $N
done

for T in 2 4 ; do
	for N in 200000 2000000 20000000 200000000 2000000000 ; do
		echo  -n N = $N T = $T ' \t '
		3/pares $T $N
	done
done
echo FIN
