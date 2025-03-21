#!/bin/bash

# Script: run_performance.sh
# Purpose: Gather performance data (timing) for the forest fire model.
#
# Part 2 of the project wants us to:
#  - Not focus on model behavior, only on time to complete.
#  - Test 3 grid sizes: N=50, N=100, N=500.
#  - Use p=0.6, M=50.
#  - Vary the number of MPI tasks to see how runtime scales.
#  - Present speedup, walltime, etc.

# Remove any old data file
rm -f performance_data.txt

# Grid sizes to test
grid_sizes="50 100 500"

# p=0.6 and M=50
p=0.6
M=50

# Range of MPI ranks to test
# Adjust as you see fit for your cluster.
mpi_tasks="1 2 4 8 16"

# Compile your code if needed (uncomment if you want automatic compile)
# mpicxx -O3 -o forestFire forestFire.cpp

echo "Starting performance runs..."
for N in $grid_sizes
do
    for np in $mpi_tasks
    do
        # Run the code with the chosen N, p, M, and number of ranks
        # forestFire prints:  nproc  N  p  M  avgSteps  avgTime  bottomFraction
        OUTPUT=$( mpirun -n $np ./forestFire $N $p $M )
        
        # We'll prepend the chosen N, p, M, and np for clarity
        # so each line in performance_data.txt looks like:
        #   N p M np  nproc  N  p  M  avgSteps  avgTime  bottomFraction
        echo "$N $p $M $np $OUTPUT" >> performance_data.txt
        
        echo "Done with N=$N, np=$np => $OUTPUT"
    done
done

echo "All runs complete. Data saved to performance_data.txt"
