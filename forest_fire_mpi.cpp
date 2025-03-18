#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>


// Define states
enum states {
    EMPTY = 0,
    TREE = 1,
    BURNING = 2,
    DEAD = 3
};

// Distribute rows over mpi tasks
void distributeRows(int N, int iproc, int nproc, int &rStart, int &rEnd) {
    int rankRows = N / nproc;
    rStart = iproc * rankRows;
    if (iproc == nproc - 1)
        rEnd = N;
    else
        rEnd = rStart + rankRows;
}

// Initialize local grid
std::vector<std::vector<states>> initLocalGrid(int globalN, double p, int rStart, int rEnd, int iproc, unsigned int seed) {
    int localRows = rEnd - rStart;

    // Add 2 halo rows
    std::vector<std::vector<states>> grid(localRows + 2, std::vector<states>(globalN, EMPTY));
    
    srand(seed);

    // Fill rows (except halos) with 'TREE' with probs p
    for (int i = 1; i <= localRows; i++) {
        for (int j = 0; j < globalN; j++) {
            double r = double(rand()) / RAND_MAX;
            if (r < p) {
                grid[i][j] = TREE;
            } else {
                grid[i][j] = EMPTY;
            }
        }
    }
    // Ignite first row of trees (if included in this rank)
    if (rStart == 0) {
        for (int j = 0; j < globalN; j++) {
            if (grid[1][j] == TREE)
                grid[1][j] = BURNING;
        }
    }
    return grid;
}

// Exchange halo rows between neighboring MPI tasks
void exchangeBoundaries(std::vector<std::vector<states>> &grid, int rStart, int rEnd, int iproc, int nproc) {
    int globalN = grid[0].size();
    int localRows = rEnd - rStart;

    // Finds rank of neighbouring process or sets NULL
    int high = (iproc == 0) ? MPI_PROC_NULL : iproc - 1;
    int low = (iproc == nproc - 1) ? MPI_PROC_NULL : iproc + 1;
    
    // Send first row to halo row 0
    MPI_Send(grid[1].data(), globalN, MPI_INT, high, 0, MPI_COMM_WORLD);
    MPI_Recv(grid[0].data(), globalN, MPI_INT, high, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // Send last row to halo row localRows+1
    MPI_Send(grid[localRows].data(), globalN, MPI_INT, low, 1, MPI_COMM_WORLD);
    MPI_Recv(grid[localRows+1].data(), globalN, MPI_INT, low, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// One time step and returns false if fire extinguished
bool step(std::vector<std::vector<states>> &grid, int rStart, int rEnd) {
    int globalN = grid[0].size();
    int localRows = rEnd - rStart;

    // Copy of current grid
    std::vector<std::vector<states>> newGrid = grid;
    bool anyBurning = false;
    
    // Loop over rows (not halso)
    for (int i = 1; i <= localRows; i++) {
        for (int j = 0; j < globalN; j++) {
            if (grid[i][j] == BURNING) {
                newGrid[i][j] = DEAD;
                
                // Spread fire to neighbours
                int di[4] = {-1, 1, 0, 0};
                int dj[4] = {0, 0, -1, 1};
                for (int k = 0; k < 4; k++) {
                    int ni = i + di[k];
                    int nj = j + dj[k];
                    if (ni >= 0 && ni <= localRows + 1 && nj >= 0 && nj < globalN) {
                        if (grid[ni][nj] == TREE)
                            newGrid[ni][nj] = BURNING;
                    }
                }
            }
        }
    }
    grid = newGrid;
    
    // Check for burning trees
    for (int i = 1; i <= localRows; i++) {
        for (int j = 0; j < globalN; j++) {
            if (grid[i][j] == BURNING) {
                anyBurning = true;
                break;
            }
        }
        if (anyBurning)
            break;
    }
    
    return anyBurning;
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int iproc, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    // Defaults
    int N = 100;
    double p = 0.5;
    int M = 50;

    // Check cmd line for args and overwrite
    if (iproc == 0) {
        if (argc >= 2) N = std::atoi(argv[1]);
        if (argc >= 3) p = std::atof(argv[2]);
        if (argc >= 4) M = std::atoi(argv[3]);
        std::cout << "N=" << N << " p=" << p << " M=" << M 
                  << "  MPI tasks:" << nproc << std::endl;
    }

    // Broadcast to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Distribute rows
    int rStart, rEnd;
    distributeRows(N, iproc, nproc, rStart, rEnd);
    int localRows = rEnd - rStart;

    // Accumulate results over all M sims
    double totalSteps = 0.0; 
    double totalTime = 0.0;
    int totalHitBottom = 0;

    // Loop over M sim runs
    for (int run = 0; run < M; run++) {
        
        // Unique seed for each process
        unsigned int seed = static_cast<unsigned int>(time(NULL)) + run * 1000 + iproc;

        // Init local grid
        auto grid = initLocalGrid(N, p, rStart, rEnd, iproc, seed);

        double startTime = MPI_Wtime();
        int steps = 0;

        // Continue until no cells burning
        while (true) {
            exchangeBoundaries(grid, rStart, rEnd, iproc, nproc);
            bool localBurning = step(grid, rStart, rEnd);
            int localFlag = localBurning ? 1 : 0;
            int globalFlag = 0;
            MPI_Allreduce(&localFlag, &globalFlag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (globalFlag > 0)
                steps++;
            else
                break;
        }
        double runTime = MPI_Wtime() - startTime;
        totalSteps += steps;
        totalTime += runTime;
        
        // Check if the fire reached the bottom (only the rank owning the bottom row does this).
        bool bottomReachedLocal = false;
        if (rEnd == N) {
            for (int j = 0; j < N; j++) {
                if (grid[localRows][j] == DEAD || grid[localRows][j] == BURNING) {
                    bottomReachedLocal = true;
                    break;
                }
            }
        }
        int localBottom = bottomReachedLocal ? 1 : 0;
        int globalBottom = 0;
        MPI_Reduce(&localBottom, &globalBottom, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (iproc == 0 && globalBottom > 0)
            totalHitBottom++;
    }
    
    // Calculate averages.
    double avgSteps = totalSteps / M;
    double avgTime = totalTime / M;
    double bottomFraction = (iproc == 0) ? (double(totalHitBottom) / M) : 0.0;
    
    // Print one standardized result line (for use in your report).
    if (iproc == 0) {
        std::cout << "RESULT: N=" << N 
                  << " p=" << p 
                  << " M=" << M 
                  << " avg_steps=" << avgSteps 
                  << " avg_time=" << avgTime 
                  << " bottom_fraction=" << bottomFraction 
                  << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}
