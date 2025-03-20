#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>

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

// Read initial global grid from a file
// File format: first integer is grid size (N) (grid is N x N), followed by N*N integers representing states.
std::vector<std::vector<states>> readInitialGrid(const std::string &filename, int &globalN) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    infile >> globalN;
    std::vector<std::vector<states>> grid(globalN, std::vector<states>(globalN, EMPTY));
    for (int i = 0; i < globalN; i++) {
        for (int j = 0; j < globalN; j++) {
            int val;
            infile >> val;
            grid[i][j] = static_cast<states>(val);
        }
    }
    infile.close();
    return grid;
}

// Initialize local grid from file data
std::vector<std::vector<states>> initLocalGridFromFile(const std::vector<std::vector<states>> &globalGrid, int rStart, int rEnd, int iproc) {
    int globalN = globalGrid.size();
    int localRows = rEnd - rStart;
    // Add 2 halo rows
    std::vector<std::vector<states>> grid(localRows + 2, std::vector<states>(globalN, EMPTY));
    // Copy the rows from the global grid into the local grid (into indices 1..localRows)
    for (int i = 0; i < localRows; i++) {
        grid[i+1] = globalGrid[rStart + i];
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

// Initialize local grid (random generation)
std::vector<std::vector<states>> initLocalGrid(int globalN, double p, int rStart, int rEnd, int iproc, unsigned int seed) {
    int localRows = rEnd - rStart;
    // Add 2 halo rows
    std::vector<std::vector<states>> grid(localRows + 2, std::vector<states>(globalN, EMPTY));
    
    srand(seed);
    // Fill rows (except halos) with 'TREE' with probability p
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
    // Loop over rows (not halos)
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
    bool useFile = false;
    std::string filename = "";
    
    // Check cmd line for args and overwrite
    // Expected usage: ./program [N] [p] [M] [optional: filename]
    if (iproc == 0) {
        if (argc >= 2) N = std::atoi(argv[1]);
        if (argc >= 3) p = std::atof(argv[2]);
        if (argc >= 4) M = std::atoi(argv[3]);
        if (argc >= 5) {
            useFile = true;
            filename = argv[4];
        }
        std::cout << "N=" << N << " p=" << p << " M=" << M 
                  << "  MPI tasks:" << nproc 
                  << (useFile ? " Using file: " + filename : " Using random grid") 
                  << std::endl;
    }
    
    // Broadcast parameters to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int useFileInt = useFile ? 1 : 0;
    MPI_Bcast(&useFileInt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    useFile = (useFileInt == 1);
    
    // If using file, broadcast filename length and filename
    int filenameLen = 0;
    if (iproc == 0) {
        filenameLen = filename.size();
    }
    MPI_Bcast(&filenameLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
    char fnameBuf[256];
    if (iproc == 0) {
        std::strcpy(fnameBuf, filename.c_str());
    }
    MPI_Bcast(fnameBuf, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (iproc != 0) {
        filename = std::string(fnameBuf);
    }
    
    // If using a file, read and broadcast the grid size (N) and the grid
    std::vector<std::vector<states>> globalGrid;
    if (useFile) {
        if (iproc == 0) {
            globalGrid = readInitialGrid(filename, N);  // N is updated here from the file
        }
        // Broadcast the updated N to all processes so that they use the file's grid size
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // Broadcast the global grid by flattening it
        std::vector<int> flatGrid;
        if (iproc == 0) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    flatGrid.push_back(static_cast<int>(globalGrid[i][j]));
                }
            }
        } else {
            flatGrid.resize(N * N);
        }
        MPI_Bcast(flatGrid.data(), N * N, MPI_INT, 0, MPI_COMM_WORLD);
        // Reconstruct globalGrid on all processes
        if (iproc != 0) {
            globalGrid.resize(N, std::vector<states>(N, EMPTY));
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    globalGrid[i][j] = static_cast<states>(flatGrid[i * N + j]);
                }
            }
        }
    }
    
    // Distribute rows
    int rStart, rEnd;
    distributeRows(N, iproc, nproc, rStart, rEnd);
    int localRows = rEnd - rStart;
    
    // Accumulate results over all M simulations
    double totalSteps = 0.0; 
    double totalTime = 0.0;
    int totalHitBottom = 0;
    
    // Loop over M simulation runs
    for (int run = 0; run < M; run++) {
        // Unique seed for each process
        unsigned int seed = static_cast<unsigned int>(time(NULL)) + run * 1000 + iproc;
        
        // Init local grid: use file reading if enabled, else random generation
        std::vector<std::vector<states>> grid;
        if (useFile) {
            grid = initLocalGridFromFile(globalGrid, rStart, rEnd, iproc);
        } else {
            grid = initLocalGrid(N, p, rStart, rEnd, iproc, seed);
        }
        
        double startTime = MPI_Wtime();
        int steps = 0;
        
        // Continue simulation until no cells are burning
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
        
        // For the process owning the bottom of the grid, check if the last row caught fire
        bool bottomHitLocal = false;
        if (rEnd == N) {
            for (int j = 0; j < N; j++) {
                if (grid[localRows][j] == DEAD || grid[localRows][j] == BURNING) {
                    bottomHitLocal = true;
                    break;
                }
            }
        }
        int localBottom = bottomHitLocal ? 1 : 0;
        int globalBottom = 0;
        MPI_Reduce(&localBottom, &globalBottom, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (iproc == 0 && globalBottom > 0)
            totalHitBottom++;
    }
    
    // Average over the M simulations
    double avgSteps = totalSteps / M;
    double avgTime = totalTime / M;
    double bottomFraction = (iproc == 0) ? (double(totalHitBottom) / M) : 0.0;
    
    // Print results
    if (iproc == 0) {
        std::cout << "For N=" << N 
                  << " p=" << p 
                  << " M=" << M 
                  << ", avg_steps=" << avgSteps 
                  << " avg_time=" << avgTime 
                  << " bottom_fraction=" << bottomFraction 
                  << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}