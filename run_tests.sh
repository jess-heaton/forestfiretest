#!/bin/bash
# run_tests.sh
# This script runs the forest fire simulation using the Moore (8-neighbor) feature.
# It tests different grid sizes (N) and tree probabilities (p).
# M (number of runs) is fixed at 50.
# Adjust the number of MPI processes (-np 4) if desired.

# Array of grid sizes to test
grid_sizes=(50 100 500)

# Array of tree probabilities to test
probs=(0.4 0.6 0.8)

# Number of simulation runs (M)
M=50

# Loop over grid sizes and probabilities
for N in "${grid_sizes[@]}"; do
  for p in "${probs[@]}"; do
    echo "=============================================="
    echo "Running simulation: N=${N}, p=${p}, M=${M} (Moore, 8-neighbor)"
    echo "----------------------------------------------"
    # Run the simulation with 4 MPI processes
    mpirun -np 4 ./forest_fire_mpi $N $p $M
    echo "=============================================="
    echo ""
  done
done