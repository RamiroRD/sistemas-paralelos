#!/bin/sh

set -e

make $1

for i in 32 64 128 256 512 1024 2048 ; do
	>&2 echo "N = $i:"
	./$1 $i | cut -f 6 -d \  | sed -n 'p;n' | sed '$!N;s/\n/ \& /'
	>&2 echo ""
done
