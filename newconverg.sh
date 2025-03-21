#!/bin/bash

# Script: run_convergence.sh
# Purpose: Gather data to study model convergence for the Forest Fire code
# Usage: 
#   1) Make sure you have compiled your newforest executable 
#      (the lecture-style code that prints one line: nproc N p M avgSteps avgTime bottomFraction).
#   2) chmod +x run_convergence.sh
#   3) ./run_convergence.sh
#
# This script will:
#  - fix N=100
#  - vary p in a small range
#  - vary M from 10 up to 100 in increments of 10
#  - run the code on 1 MPI rank
#  - append results to convergence_data.txt

# Remove any old data file
rm -f convergence_data.txt

# Define the probabilities you want to test
p_values="0.3 0.4 0.5 0.6 0.7"

# Define the M values for repeat runs
M_values="10 20 30 40 50 60 70 80 90 100"

# Choose how many MPI ranks to run on (1 is simplest for a pure convergence check)
NP=1

# Loop over probabilities and M values
for p in $p_values
do
    for M in $M_values
    do
        # Run your newforest code, capturing its single-line output
        # Syntax: mpirun -n $NP ./newforest <N> <p> <M>
        OUTPUT=$( mpirun -n $NP ./newforest 100 $p $M )

        # The code prints something like:  "1 100 0.4 30 25.7 0.74 0.62"
        # We'll prepend 'p' and 'M' so we have them explicitly, 
        # in case you want them in your output file:
        # e.g. "0.4 30 1 100 0.4 30 25.7 0.74 0.62"

        echo "$p $M $OUTPUT" >> convergence_data.txt

        # You can echo to the screen for progress:
        echo "Done with p=$p, M=$M => $OUTPUT"
    done
done

echo "All runs complete. Data saved to convergence_data.txt"
