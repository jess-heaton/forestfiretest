#!/bin/bash
# run_convergence.sh
# This script runs convergence tests for a fixed grid size N = 100.
# It varies the tree probability p and the number of independent runs M.
# The results (avg_steps, avg_time, and bottom_fraction) are written to convergence_results.txt.

output_file="convergence_results.txt"
echo "p M avg_steps avg_time bottom_fraction" > $output_file

# Loop over different p values.
for p in 0.2 0.4 0.6 0.8; do
    # Loop over different values of M.
    for M in 10 25 50 75 100; do
        echo "Running for p=$p and M=$M"
        # Run the simulation with fixed N=100 and, for example, 4 MPI processes.
        mpirun -np 4 ./forest_fire_mpi 100 $p $M >> $output_file
    done
done

echo "Convergence test results written to $output_file"