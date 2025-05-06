from mpi4py import MPI
import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
import seaborn as sns
import time
import pandas as pd

# MPI Initialization
comm = MPI.COMM_WORLD
rank = comm.Get_rank()  # Process rank
size = comm.Get_size()  # Total number of processes

# Constants
kB = 1.3806504e-23  # Boltzmann constant [J/K]
eV2J = 1.60218e-19  # eV to Joules
A2m = 1.0e-10  # Angstroms to meters
ps2s = 1.0e-12  # Picoseconds to seconds

# Conversion factor
convert = (eV2J**2) / (ps2s * A2m)

# Simulation parameters
T = 500  # Temperature [K]
dt = 0.000766  # Timestep [ps]
V = 109**3  # Volume in cubic angstroms

# Define sample intervals and number of samples for testing
s_values = np.arange(50, 2000, 50)
p_values = np.arange(50, 2000, 50)

# Master node loads data and distributes work
if rank == 0:
    flux_data = np.loadtxt("heat_flux.dat")
    timesteps = flux_data[:, 0]
    Jx = flux_data[:, 1]
    Jy = flux_data[:, 2]
    Jz = flux_data[:, 3]

    # Define truncation time in picoseconds
    truncation_time = 1000
    truncation_timestep = int(truncation_time / dt)

    # Print truncation details
    print(f"Truncation time: {truncation_time} ps")
    print(f"Truncation timestep: {truncation_timestep} steps")
    print(f"Timestep, dt: {dt} ps")

    # Truncate data
    timesteps = timesteps[timesteps <= truncation_timestep]
    Jx = Jx[:len(timesteps)]
    Jy = Jy[:len(timesteps)]
    Jz = Jz[:len(timesteps)]

    # Create list of (s, p) combinations
    task_list = [(s, p) for s in s_values for p in p_values]

    # Distribute tasks to workers
    task_chunks = np.array_split(task_list, size)
else:
    flux_data = None
    timesteps = None
    Jx = None
    Jy = None
    Jz = None
    task_chunks = None

# Broadcast data to all workers
timesteps = comm.bcast(timesteps, root=0)
Jx = comm.bcast(Jx, root=0)
Jy = comm.bcast(Jy, root=0)
Jz = comm.bcast(Jz, root=0)

# Scatter task chunks to each worker
my_tasks = comm.scatter(task_chunks, root=0)

# Function to compute autocorrelation function (ACF)
def compute_acf(J, p, s, d_timestep):
    C_vals = np.zeros(p)
    counts = np.zeros(p)

    n_datapoints = int(d_timestep / s)
    if n_datapoints != len(J):
        raise RuntimeError("Calculation must use the final timestep.")

    for p_num in range(0, p):
        lag = p_num * s
        lag_downsampled = p_num
        n_windows = n_datapoints - p_num - p + 1

        if n_windows > 0:
            C_vals[p_num] = 0.0
            for i in range(n_windows):
                loc0 = i
                loct = loc0 + lag_downsampled
                J_0 = J[loc0:loc0 + p]
                J_t = J[loct:loct + p]
                C_vals[p_num] += (1 / p) * np.dot(J_0, J_t)
            C_vals[p_num] /= n_windows
            counts[p_num] = n_windows
    return C_vals, counts

# Worker nodes process assigned (s, p) pairs
thermal_data = []

for s, p in my_tasks:
    d = p * s

    timesteps_downsampled = timesteps[::s]
    Jx_downsampled = Jx[::s]
    Jy_downsampled = Jy[::s]
    Jz_downsampled = Jz[::s]

    if (int(d / s) < len(timesteps_downsampled)):
        last_d_timestep = len(timesteps_downsampled) * s

        C_x, _ = compute_acf(Jx_downsampled, p, s, last_d_timestep)
        C_y, _ = compute_acf(Jy_downsampled, p, s, last_d_timestep)
        C_z, _ = compute_acf(Jz_downsampled, p, s, last_d_timestep)

        scale = (convert / (kB * T * T * V)) * s * dt
        kappa_x = trapezoid(C_x) * scale
        kappa_y = trapezoid(C_y) * scale
        kappa_z = trapezoid(C_z) * scale


        kappa_iso = (kappa_x + kappa_y + kappa_z) / 3.0
        thermal_data.append([s, p, kappa_x, kappa_y, kappa_z])

# Gather results from all workers
all_thermal_data = comm.gather(thermal_data, root=0)

# Master node processes results
if rank == 0:
    # Flatten list of lists
    final_data = [item for sublist in all_thermal_data for item in sublist]

    # Save computed data
    np.savetxt("thermal_conductivity_data.dat", final_data, fmt="%.6f", header="s p kappa_x kappa_y kappa_z")

    # Load the computed data
    thermal_data = np.loadtxt("thermal_conductivity_data.dat")
    s_values = thermal_data[:, 0]
    p_values = thermal_data[:, 1]
    kappa_x = thermal_data[:, 2]
    kappa_y = thermal_data[:, 3]
    kappa_z = thermal_data[:, 4]

    # Compute average thermal conductivity
    kappa_avg = (kappa_x + kappa_y + kappa_z) / 3.0

    # Create heatmap grid
    unique_s = np.unique(s_values)
    unique_p = np.unique(p_values)
    kappa_matrix = np.zeros((len(unique_p), len(unique_s)))

    for i, p in enumerate(unique_p):
        for j, s in enumerate(unique_s):
            mask = (s_values == s) & (p_values == p)
            if np.any(mask):
                kappa_matrix[i, j] = np.mean(kappa_avg[mask])
                
                
    # Create a DataFrame from the kappa_matrix
    df_kappa = pd.DataFrame(kappa_matrix, index=unique_p, columns=unique_s)

    # Save to CSV
    df_kappa.to_csv("kappa_matrix.csv", index=True, header=True)

    # Plot heatmap
    plt.figure(figsize=(10, 6))
    show_values = True  # ← Toggle this to False if you don't want numbers inside the squares

    sns.heatmap(
        kappa_matrix,
        xticklabels=unique_s,
        yticklabels=unique_p,
        cmap="coolwarm",
        annot=show_values,
        fmt=".2f",             # Format numbers to 2 decimal places
        annot_kws={"size": 6}  # Smaller font to prevent overflow
    )

    plt.xlabel("Sample Interval (s)")
    plt.ylabel("Number of Samples (p)")
    plt.title("Thermal Conductivity calculated by Green-Kubo Method")
    plt.tight_layout()
    plt.savefig("kappa_heatmap.png")


    print("Finished Parallel Computation!")
