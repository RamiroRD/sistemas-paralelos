#!/bin/sh

set -e

cc -O0 -o multBloques multBloques.c

for r in 32 64 128 256 512 1024; do
	for n in 1 2 3 4; do
		>&2 echo "r = $r	n=$n"
		./multBloques $n $r 0
	done
done

