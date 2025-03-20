import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


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
s_values = np.arange(10, 110, 10)  # Example intervals
p_values = np.arange(10, 110, 10)  # Example number of samples





# Load flux data
flux_data = np.loadtxt("heat_flux.dat")
timesteps = flux_data[:, 0]
Jx = flux_data[:, 1]
Jy = flux_data[:, 2]
Jz = flux_data[:, 3]

# TRUNCATE DATA IF NECESSARY (e.g if run for too long and kappa begins drifting)
# Define truncation time in picoseconds if slicing heat_flux
truncation_time = 1000 

# Convert truncation time to timesteps
truncation_timestep = int(truncation_time / dt)

# Print truncation details
print(f"Truncation time: {truncation_time} ps")
print(f"Truncation timestep: {truncation_timestep} steps")
print(f"Timestep, dt: {dt} ps\n")

# Truncate data
timesteps = timesteps[timesteps <= truncation_timestep]
Jx = Jx[:len(timesteps)]
Jy = Jy[:len(timesteps)]
Jz = Jz[:len(timesteps)]

# Multiply by V if already scaled by volume (code is for non-scaled fluxes)
# Jx = V * Jx
# Jy = V * Jy
# Jz = V * Jz







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
        
        if (n_datapoints != len(J)):
            raise RuntimeError("Calculation of thermal conductivity must be with FINAL timestep")
        
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
        
        timesteps_downsampled = timesteps[::s]
        Jx_downsampled = Jx[::s]
        Jy_downsampled = Jy[::s]
        Jz_downsampled = Jz[::s]  

        
        if (int(d/s) < len(timesteps_downsampled)):
        
            last_d_timestep = len(timesteps_downsampled) * s
        
            C_x, _ = compute_acf(Jx_downsampled, p, s, last_d_timestep)
            C_y, _ = compute_acf(Jy_downsampled, p, s, last_d_timestep)
            C_z, _ = compute_acf(Jz_downsampled, p, s, last_d_timestep)
            
            scale = (convert / (kB * T * T * V)) * s * dt
            kappa_x = np.trapz(C_x) * scale
            kappa_y = np.trapz(C_y) * scale
            kappa_z = np.trapz(C_z) * scale
            
            thermal_data.append([s, p, kappa_x, kappa_y, kappa_z])
            
            kappa_iso = (kappa_x + kappa_y + kappa_z) / 3.0
            print(f"Thermal conductivity computed for s = {s}, p = {p}. kappa ={kappa_iso:.2f} W/mK\n")





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



# Create a DataFrame from the kappa_matrix
df_kappa = pd.DataFrame(kappa_matrix, index=unique_p, columns=unique_s)

# Save to CSV
df_kappa.to_csv("kappa_matrix.csv", index=True, header=True)


# Plot heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(kappa_matrix, xticklabels=unique_s, yticklabels=unique_p, cmap="coolwarm")
plt.xlabel("Sample Interval (s)")
plt.ylabel("Number of Samples (p)")
plt.title("Thermal Conductivity calculated by Green-Kubo Method")
plt.tight_layout()
plt.savefig("kappa_heatmap.png")
