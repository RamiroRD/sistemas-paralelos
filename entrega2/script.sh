#!/bin/bash
make

set -o xtrace

# Secuencial
for size in 512 1024 2048
do
        ./secuencial $size
done

# MPI 1 PC, 4 threads
for size in 512 1024 2048
do
        mpirun -np 4 mpi $size
done

# MPI 2 PC, 4 threads (2 c/u)
for size in 512 1024 2048
do
        mpirun -np 4 -machinefile maquinas_2_4 mpi $size
done

# MPI 2 PC, 8 threads (4 c/u)
for size in 512 1024 2048
do
        mpirun -np 8 -machinefile maquinas_2_8 mpi $size
done

# MPI-OMP 1 PC, 1 slot, 4 threads/slot
for size in 512 1024 2048
do
        mpirun -np 1 --bind-to none hibrido $size 4
done

# MPI-OMP 2 PC, 2 slots, 2 threads/slot
for size in 512 1024 2048
do
        mpirun -np 2 -machinefile maquinas_2_2 --bind-to none hibrido $size 2
done

# MPI-OMP 2 PC, 2 slots, 4 threads/slot
for size in 512 1024 2048
do
        mpirun -np 2 -machinefile maquinas_2_2 --bind-to none hibrido $size 4
done
