import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load flux data
flux_data = np.loadtxt("heat_flux.dat")
timesteps = flux_data[:, 0]
Jx = flux_data[:, 1]
Jy = flux_data[:, 2]
Jz = flux_data[:, 3]

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
s_values = np.arange(25, 550, 25)  # Example intervals
p_values = np.arange(100, 2100, 100)  # Example number of samples

# Compute ACF function
def compute_acf(J, p, s, d_timestep):

    # Compute the autocorrelation function C_J(Δt) at every d timestep.

    # J: The downsampled heat flux array.
    # p: Number of correlation samples.
    # s: Sample interval.
    # d_timestep: Current time step, must be a multiple of (p * s)
    # C_vals: numpy array of current C(0), C(s), C(2s), ..., C(d - s). len(C_vals) should equal p.
    # counts: numpy array of number of windows contributing to each corresponding C, is also of size p.
    
    C_vals = np.zeros(p)
    counts = np.zeros(p)

    for p_num in range(0, p):
    
        lag = p_num * s
        lag_downsampled = p_num

        # Find number of windows available for the current C(lag)
        n_datapoints = int(d_timestep/s)
        
        n_windows = n_datapoints - p_num - p + 1
        

        if n_windows > 0:
            C_vals[p_num] = 0.0

            for i in range(n_windows):
                loc0 = i
                loct = loc0 + lag_downsampled

                J_0 = J[loc0:loc0 + p]
                J_t = J[loct:loct + p]

                C_vals[p_num] += (1/p) * np.dot(J_0, J_t)

            C_vals[p_num] /= n_windows  # Average over available windows         
            counts[p_num] = n_windows  # Track number of contributions

    return C_vals, counts

# Compute and store thermal conductivity data
thermal_data = []

for s in s_values:
    for p in p_values:
        d = p * s
        total_timesteps = len(timesteps)
        
        if d < total_timesteps:
            C_x, _ = compute_acf(Jx, p, s, total_timesteps)
            C_y, _ = compute_acf(Jy, p, s, total_timesteps)
            C_z, _ = compute_acf(Jz, p, s, total_timesteps)
            
            scale = (convert / (kB * T * T * V)) * s * dt
            kappa_x = np.trapz(C_x) * scale
            kappa_y = np.trapz(C_y) * scale
            kappa_z = np.trapz(C_z) * scale
            
            thermal_data.append([s, p, kappa_x, kappa_y, kappa_z])

# Save computed data
np.savetxt("thermal_conductivity_data.dat", thermal_data, fmt="%.6f", header="s p kappa_x kappa_y kappa_z")

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

# Plot heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(kappa_matrix, xticklabels=unique_s, yticklabels=unique_p, cmap="coolwarm", annot=True, fmt=".2f")
plt.xlabel("Sample Interval (s)")
plt.ylabel("Number of Samples (p)")
plt.title("Heat Map of Thermal Conductivity")
plt.savefig("kappa_heatmap.png")
