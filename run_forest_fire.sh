#!/bin/bash
#SBATCH --job-name=forest_fire_perf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
#SBATCH --partition=teach_cpu     # or another valid partition
#SBATCH --output=forest_fire_%j.out

module load mpi

p=0.6
M=50

for N in 50 100 500; do
    echo "Running with grid size N=${N}"
    mpirun -np 8 ./your_executable ${N} ${p} ${M}
done
