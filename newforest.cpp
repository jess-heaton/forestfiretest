#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>        // C++11 random
#include <string>
#include <cstdlib>       // For std::atoi, std::atof
#include <ctime>         // For time(NULL)
#include <algorithm>     // For std::max, std::min
#include <sstream>       // For std::stringstream

/*****************************************************************************
 *                        1. STATE DEFINITIONS
 *****************************************************************************/
static const int EMPTY   = 0;
static const int TREE    = 1;
static const int BURNING = 2;
static const int DEAD    = 3;

/*****************************************************************************
 *                     2. DISTRIBUTE GRID (LECTURE STYLE)
 *****************************************************************************/
void distribute_grid(int N, int iproc, int nproc, int &i0, int &i1)
{
    // Simple block decomposition:
    int blockSize = N / nproc;
    i0 = iproc * blockSize;
    i1 = (iproc == nproc - 1) ? N : i0 + blockSize;
}

/*****************************************************************************
 *  3. READING THE GRID FROM FILE *WITHOUT* N AT THE TOP,
 *     THEN BROADCASTING (LECTURE STYLE)
 *
 *  We read all lines, parse them into an NxN flattened array, and deduce N
 *  automatically (assuming it's square).
 *****************************************************************************/
std::vector<int> read_and_broadcast_grid_noSize(const std::string &filename,
                                                int &N, int iproc, int nproc)
{
    std::vector<int> flatData; // Will store the entire grid flattened.

    if (iproc == 0) {
        // Rank 0 reads
        std::ifstream infile(filename);
        if (!infile) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // We'll read all lines, parse them as integers, store in flatData
        std::string line;
        std::vector<std::vector<int>> rows;
        while (std::getline(infile, line)) {
            if (line.empty()) continue; // skip blank lines if any
            std::stringstream ss(line);
            std::vector<int> rowVals;
            int val;
            while (ss >> val) {
                rowVals.push_back(val);
            }
            if (!rowVals.empty()) {
                rows.push_back(rowVals);
            }
        }
        infile.close();

        // Check rowCount, colCount
        int rowCount = (int)rows.size();
        if (rowCount == 0) {
            std::cerr << "Error: The file is empty or invalid: " << filename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        int colCount = (int)rows[0].size();
        for (int r = 1; r < rowCount; r++) {
            if ((int)rows[r].size() != colCount) {
                std::cerr << "Error: Inconsistent number of columns in " << filename << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        if (rowCount != colCount) {
            std::cerr << "Error: The file does not describe a square grid ("
                      << rowCount << "x" << colCount << ")" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        N = rowCount; // deduce N

        // Flatten row-by-row
        flatData.resize(N * N);
        for (int r = 0; r < N; r++) {
            for (int c = 0; c < N; c++) {
                flatData[r*N + c] = rows[r][c];
            }
        }
    }

    // First broadcast the deduced N
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now broadcast the entire NxN data
    if (iproc != 0) {
        flatData.resize(N * N);
    }
    MPI_Bcast(flatData.data(), N*N, MPI_INT, 0, MPI_COMM_WORLD);

    return flatData;
}

/*****************************************************************************
 *          4. RANDOM GRID INITIALIZATION (LECTURE-STYLE W/ <random>)
 *****************************************************************************/
std::vector<int> generate_random_grid(int N, double p, int iproc, int nproc,
                                      int run, unsigned int baseSeed=42)
{
    // We'll only generate the random grid on rank 0, then broadcast it
    std::vector<int> gridData(N * N, EMPTY);

    if (iproc == 0) {
        unsigned int seed = baseSeed + 1000 * run;
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (int i = 0; i < N*N; i++) {
            double r = dist(gen);
            if (r < p) {
                gridData[i] = TREE;
            } else {
                gridData[i] = EMPTY;
            }
        }
        // Ignite top row
        for (int j = 0; j < N; j++) {
            if (gridData[j] == TREE) {
                gridData[j] = BURNING;
            }
        }
    }

    MPI_Bcast(gridData.data(), N*N, MPI_INT, 0, MPI_COMM_WORLD);
    return gridData;
}

/*****************************************************************************
 *          5. STEP FUNCTION (LECTURE-STYLE W/ PARTIAL UPDATES + ALLREDUCE)
 *****************************************************************************/
bool step_forest_fire(std::vector<int> &oldGridData,
                      int N, int iproc, int nproc,
                      int i0, int i1, bool moore=false)
{
    std::vector<int> newGridData(oldGridData); 
    bool localBurning = false;

    // Offsets for neighbors
    int d4i[4] = {-1,  1,  0,  0};
    int d4j[4] = { 0,  0, -1,  1};

    int d8i[8] = {-1, -1, -1,  0,  0,  1,  1,  1};
    int d8j[8] = {-1,  0,  1, -1,  1, -1,  0,  1};

    // Process only rows [i0..i1)
    for (int i = i0; i < i1; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i*N + j;
            if (oldGridData[idx] == BURNING) {
                // BURNING -> DEAD
                newGridData[idx] = DEAD;

                int nCount = (moore ? 8 : 4);
                for (int k = 0; k < nCount; k++) {
                    int ni = i + (moore ? d8i[k] : d4i[k]);
                    int nj = j + (moore ? d8j[k] : d4j[k]);
                    if (ni >= 0 && ni < N && nj >= 0 && nj < N) {
                        int nidx = ni*N + nj;
                        if (oldGridData[nidx] == TREE) {
                            newGridData[nidx] = BURNING;
                        }
                    }
                }
            }
        }
    }

    // Copy partial updates into oldGridData (local slice)
    for (int i = i0; i < i1; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i*N + j;
            oldGridData[idx] = newGridData[idx];
        }
    }

    // Combine partial updates across all ranks
    MPI_Allreduce(MPI_IN_PLACE, oldGridData.data(), N*N, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // Check if any cells are still BURNING
    for (int i = i0; i < i1; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i*N + j;
            if (oldGridData[idx] == BURNING) {
                localBurning = true;
                break;
            }
        }
        if (localBurning) break;
    }

    int localFlag = localBurning ? 1 : 0;
    int globalFlag;
    MPI_Allreduce(&localFlag, &globalFlag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return (globalFlag > 0);
}

/*****************************************************************************
 *                  6. MAIN PROGRAM
 *****************************************************************************/
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int iproc, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Default parameters
    int    N       = 100;
    double p       = 0.5;
    int    M       = 50;
    bool   useFile = false;
    std::string filename;

    // Parse command line (rank 0 only)
    if (iproc == 0) {
        // e.g. ./program 100 0.5 50 myGrid.txt
        if (argc >= 2) N = std::atoi(argv[1]);
        if (argc >= 3) p = std::atof(argv[2]);
        if (argc >= 4) M = std::atoi(argv[3]);
        if (argc >= 5) {
            useFile = true;
            filename = argv[4];
        }
    }

    // Broadcast these
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int useFileInt = (useFile ? 1 : 0);
    MPI_Bcast(&useFileInt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    useFile = (useFileInt == 1);

    // Distribute domain among ranks
    int i0, i1;
    distribute_grid(N, iproc, nproc, i0, i1);

    // We'll accumulate results over M runs
    double totalSteps = 0.0;
    double totalTime  = 0.0;
    int    totalHitBottom = 0;

    // If reading from file, read entire NxN *without* leading size, broadcast to all.
    // Then store in fileGridData for re-use each run.
    std::vector<int> fileGridData;
    if (useFile) {
        if (iproc == 0 && !filename.empty()) {
            std::cerr << "Reading grid from " << filename 
                      << " (no size on first line, deducing N...)\n";
        }
        fileGridData = read_and_broadcast_grid_noSize(filename, N, iproc, nproc);

        // re-distribute domain in case N changed
        distribute_grid(N, iproc, nproc, i0, i1);
    }

    // Main loop over M simulations
    for (int run = 0; run < M; run++)
    {
        // Build initial grid
        std::vector<int> gridData;
        if (useFile) {
            // We have fileGridData => copy it
            gridData = fileGridData;

            // Also ignite top row if it's not already burning
            if (iproc == 0) {
                for (int j = 0; j < N; j++) {
                    int idx = j;
                    if (gridData[idx] == TREE) {
                        gridData[idx] = BURNING;
                    }
                }
            }
            // Make sure all ranks see that top row
            MPI_Bcast(gridData.data(), N*N, MPI_INT, 0, MPI_COMM_WORLD);
        } else {
            // Random generation
            gridData = generate_random_grid(N, p, iproc, nproc, run);
        }

        double startTime = MPI_Wtime();
        int steps = 0;

        // Keep stepping until no cells are burning
        while (true) {
            bool burning = step_forest_fire(gridData, N, iproc, nproc, i0, i1);
            if (burning) steps++;
            else break;
        }

        double runTime = MPI_Wtime() - startTime;
        totalSteps += steps;
        totalTime  += runTime;

        // Check if bottom row is burnt or burning
        bool bottomHitLocal = false;
        if (iproc == 0) {
            for (int j = 0; j < N; j++) {
                int idx = (N - 1) * N + j;
                if (gridData[idx] == DEAD || gridData[idx] == BURNING) {
                    bottomHitLocal = true;
                    break;
                }
            }
        }
        int localFlag = bottomHitLocal ? 1 : 0;
        int globalFlag;
        MPI_Allreduce(&localFlag, &globalFlag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (iproc == 0 && globalFlag > 0) {
            totalHitBottom++;
        }
    }

    // Compute final averages
    double avgSteps = totalSteps / M;
    double avgTime  = totalTime  / M;
    double bottomFraction = 0.0;
    if (iproc == 0) {
        bottomFraction = double(totalHitBottom) / double(M);
    }

    // Print a single line of data for easy Python parsing:
    // Format: nproc N p M avgSteps avgTime bottomFraction
    if (iproc == 0) {
        std::cout << nproc << " "
                  << N << " "
                  << p << " "
                  << M << " "
                  << avgSteps << " "
                  << avgTime << " "
                  << bottomFraction
                  << std::endl;
    }

    MPI_Finalize();
    return 0;
}
