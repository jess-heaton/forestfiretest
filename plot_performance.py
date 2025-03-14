import numpy as np
import matplotlib.pyplot as plt

# Load data from the performance file.
data = np.loadtxt('performance_scaling.txt', skiprows=1)
# Columns: 0: N, 1: np, 2: p, 3: M, 4: avg_steps, 5: avg_time, 6: bottom_fraction

# Get unique grid sizes.
Ns = np.unique(data[:,0])

plt.figure(figsize=(10,6))
for N in Ns:
    subset = data[data[:,0]==N]
    plt.plot(subset[:,1], subset[:,5], marker='o', label=f'N={int(N)}')
plt.xlabel('Number of MPI Processes')
plt.ylabel('Average Runtime (s)')
plt.title('Scaling of Forest Fire Simulation')
plt.legend()
plt.grid(True)
plt.show()
