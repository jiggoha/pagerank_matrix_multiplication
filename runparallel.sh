#!/usr/bin/bash
#PBS -l procs=2,tpn=2,mem=34gb
#PBS -l walltime=1:00:00
#PBS -N Amelia
#PBS -r n
#PBS -j oe
#PBS -q cpsc424

# Load necessary module files
module load Langs/Intel/14 MPI/OpenMPI/1.6.5

# run the client for each different thread size
for n in 512 1024 2048 8192; do
	for p in 1 5 15; do
		mpiexec -n 2 ~/pagerank_matrix_multiplication/matrix_multiplication $n $p
done


#PBS -l procs=8,tpn=2,mem=136gb
#PBS -l walltime=1:00:00
#PBS -N Jay
#PBS -r n
#PBS -j oe
#PBS -q cpsc424

# Load necessary module files
module load Langs/Intel/14 MPI/OpenMPI/1.6.5

# Load necessary module files
module load Langs/Intel/14 MPI/OpenMPI/1.6.5

# run the client for each different thread size
for n in 512 1024 2048 8192; do
	for p in 1 5 15; do
		mpiexec -n 8 ~/pagerank_matrix_multiplication/matrix_multiplication $n $p
done
