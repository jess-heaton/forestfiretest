#include <iostream>
#include <vector>
#include <cstdlib>   // For rand(), RAND_MAX
#include <ctime>     // For time(nullptr)

// Possible cell states
enum CellState { EMPTY, TREE, BURNING, DEAD };

int main()
{
    // --- 1. Read user input ---
    int N;
    double p;
    std::cout << "Enter grid size N: ";
    std::cin >> N;
    std::cout << "Enter tree probability p (0 to 1): ";
    std::cin >> p;

    // --- 2. Initialize the grid ---
    // We'll store the grid as a 2D vector of CellStates.
    std::srand(static_cast<unsigned int>(std::time(nullptr))); 
    std::vector<std::vector<CellState>> grid(N, std::vector<CellState>(N, EMPTY));

    // Fill grid: each cell has probability p of being a living TREE.
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            double r = static_cast<double>(std::rand()) / RAND_MAX;
            if(r < p){
                grid[i][j] = TREE;
            } else {
                grid[i][j] = EMPTY;
            }
        }
    }

    // --- 3. Ignite the top row ---
    for(int j = 0; j < N; j++){
        if(grid[0][j] == TREE){
            grid[0][j] = BURNING;
        }
    }

    // --- 4. Simulate until no cells are BURNING ---
    int steps = 0;
    bool fireStillBurning = true;

    while(fireStillBurning){
        fireStillBurning = false; 
        // We use a copy to hold updated states for the new step
        std::vector<std::vector<CellState>> newGrid = grid;

        // For each cell that is currently BURNING, spread fire to neighbors
        // and turn it to DEAD for the next step
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                if(grid[i][j] == BURNING){
                    // This burning tree becomes DEAD
                    newGrid[i][j] = DEAD;

                    // Spread fire to up/down/left/right neighbors if they are TREE
                    int di[4] = {-1, 1, 0, 0};
                    int dj[4] = {0, 0, -1, 1};
                    for(int k = 0; k < 4; k++){
                        int ni = i + di[k];
                        int nj = j + dj[k];
                        // Check boundaries
                        if(ni >= 0 && ni < N && nj >= 0 && nj < N){
                            if(grid[ni][nj] == TREE){
                                newGrid[ni][nj] = BURNING;
                            }
                        }
                    }
                }
            }
        }

        // Check if we ignited any new cells in newGrid
        // If yes, then the fire is still burning
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                if(newGrid[i][j] == BURNING){
                    fireStillBurning = true;
                    break;
                }
            }
            if(fireStillBurning) break;
        }

        // Update grid
        grid = newGrid;
        if(fireStillBurning){
            steps++;
        }
    }

    // --- 5. Check whether fire reached bottom row ---
    bool reachedBottom = false;
    for(int j = 0; j < N; j++){
        if(grid[N - 1][j] == DEAD || grid[N - 1][j] == BURNING){
            reachedBottom = true;
            break;
        }
    }

    // --- 6. Print results ---
    std::cout << "Fire took " << steps << " steps to finish.\n";
    if(reachedBottom) {
        std::cout << "Fire *did* reach the bottom row.\n";
    } else {
        std::cout << "Fire did *not* reach the bottom row.\n";
    }

    return 0;
}
