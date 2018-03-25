#!/bin/sh

make matrices
make matrices.original

echo "Programa original:"

for i in 32 64 128 256 512 1024 ; do
	echo "matrices.original $i"
	./matrices.original $i
done

echo "Programa optimizado:"
for i in 32 64 128 256 512 1024 ; do
	echo "matrices $i"
	./matrices $i
done

