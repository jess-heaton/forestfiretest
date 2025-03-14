#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

// Define cell states
enum CellState {
    EMPTY = 0,
    TREE = 1,
    BURNING = 2,
    DEAD = 3
};

// Distribute rows among MPI tasks (1D domain decomposition)
void distributeRows(int N, int iproc, int nproc, int &rowStart, int &rowEnd) {
    int rowsPerRank = N / nproc;
    rowStart = iproc * rowsPerRank;
    if (iproc == nproc - 1)
        rowEnd = N;
    else
        rowEnd = rowStart + rowsPerRank;
}

// Initialize the local grid for one simulation run.
// The grid has extra halo rows: indices 0 and localRows+1.
std::vector<std::vector<CellState>> initLocalGrid(int globalN, double p, int rowStart, int rowEnd, int iproc, unsigned int seed) {
    int localRows = rowEnd - rowStart;
    std::vector<std::vector<CellState>> grid(localRows + 2, std::vector<CellState>(globalN, EMPTY));
    
    srand(seed); // Set process-specific random seed

    // Fill real rows (indices 1 to localRows) with TREE with probability p.
    for (int i = 1; i <= localRows; i++) {
        for (int j = 0; j < globalN; j++) {
            double r = double(rand()) / RAND_MAX;
            grid[i][j] = (r < p) ? TREE : EMPTY;
        }
    }
    // If this rank contains the global top row, ignite trees there.
    if (rowStart == 0) {
        for (int j = 0; j < globalN; j++) {
            if (grid[1][j] == TREE)
                grid[1][j] = BURNING;
        }
    }
    return grid;
}

// Exchange halo rows between neighboring MPI tasks.
void exchangeBoundaries(std::vector<std::vector<CellState>> &grid, int rowStart, int rowEnd, int iproc, int nproc) {
    int globalN = grid[0].size();
    int localRows = rowEnd - rowStart;
    int above = (iproc == 0) ? MPI_PROC_NULL : iproc - 1;
    int below = (iproc == nproc - 1) ? MPI_PROC_NULL : iproc + 1;
    
    // Send first real row upward; receive into halo row 0.
    MPI_Send(grid[1].data(), globalN, MPI_INT, above, 0, MPI_COMM_WORLD);
    MPI_Recv(grid[0].data(), globalN, MPI_INT, above, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // Send last real row downward; receive into halo row localRows+1.
    MPI_Send(grid[localRows].data(), globalN, MPI_INT, below, 1, MPI_COMM_WORLD);
    MPI_Recv(grid[localRows+1].data(), globalN, MPI_INT, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// Perform one simulation step (fire spread) and return true if any burning cells remain.
bool doOneStep(std::vector<std::vector<CellState>> &grid, int rowStart, int rowEnd) {
    int globalN = grid[0].size();
    int localRows = rowEnd - rowStart;
    std::vector<std::vector<CellState>> newGrid = grid;
    bool anyBurning = false;
    
    // Loop over real rows only.
    for (int i = 1; i <= localRows; i++) {
        for (int j = 0; j < globalN; j++) {
            if (grid[i][j] == BURNING) {
                newGrid[i][j] = DEAD; // Burning cell becomes dead.
                // Spread fire to neighbors (Von Neumann: up, down, left, right).
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
    // Check if any burning cell remains.
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
    
    // Command-line arguments:
    // argv[1] = N (grid size), argv[2] = p (tree probability), argv[3] = M (number of runs)
    int N = 100;
    double p = 0.5;
    int M = 50;
    if (iproc == 0) {
        if (argc >= 2) N = std::atoi(argv[1]);
        if (argc >= 3) p = std::atof(argv[2]);
        if (argc >= 4) M = std::atoi(argv[3]);
        std::cout << "Running Forest Fire with N=" << N << ", p=" << p << ", M=" << M 
                  << ", using " << nproc << " MPI tasks." << std::endl;
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Distribute rows among tasks.
    int rowStart, rowEnd;
    distributeRows(N, iproc, nproc, rowStart, rowEnd);
    int localRows = rowEnd - rowStart;
    
    // Variables for accumulating results over M runs.
    double totalSteps = 0.0;
    double totalTime = 0.0;
    int totalBottom = 0; // Count how many runs reached the bottom (only computed on the rank owning the bottom).
    
    for (int run = 0; run < M; run++) {
        // Use a unique seed per run and process.
        unsigned int seed = static_cast<unsigned int>(time(NULL)) + run * 1000 + iproc;
        auto grid = initLocalGrid(N, p, rowStart, rowEnd, iproc, seed);
        
        double startTime = MPI_Wtime();
        int steps = 0;
        while (true) {
            exchangeBoundaries(grid, rowStart, rowEnd, iproc, nproc);
            bool localBurning = doOneStep(grid, rowStart, rowEnd);
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
        if (rowEnd == N) {
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
            totalBottom++;
    }
    
    // Calculate averages.
    double avgSteps = totalSteps / M;
    double avgTime = totalTime / M;
    double bottomFraction = (iproc == 0) ? (double(totalBottom) / M) : 0.0;
    
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