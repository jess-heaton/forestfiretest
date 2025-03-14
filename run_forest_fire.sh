#!/bin/bash
# Example job script for BlueCrystal (adjust scheduler directives as needed)
#SBATCH --job-name=forest_fire_perf
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:30:00
#SBATCH --output=forest_fire_%j.out

# Load the MPI module (if needed)
module load mpi

# Define parameters
p=0.6
M=50

# Loop over different grid sizes
for N in 50 100 500; do
    echo "Running with grid size N=${N}"
    
    # Example: run with 8 MPI tasks
    mpirun -np 8 ./your_executable ${N} ${p} ${M}
    
    # You can add more runs with different MPI task counts for the same grid size
    # e.g., 1, 4, 8, 16, etc.
done
