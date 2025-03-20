import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import pdb


# Constants (same as in LAMMPS script)
kB = 1.3806504e-23  # Boltzmann constant [J/K]
eV2J = 1.60218e-19  # eV to Joules
A2m = 1.0e-10  # Angstroms to meters
ps2s = 1.0e-12  # Picoseconds to seconds

# Conversion factor (from LAMMPS script)
convert = (eV2J**2) / (ps2s * A2m)

# Simulation parameters (should match LAMMPS)
T = 500  # Temperature [K]
dt = 0.000766  # Timestep [ps]
p = 100  # Number of correlation samples
s = 20  # Sample interval
d = p * s  # Correlation length in number of timesteps

V = 109**3  # Volume in cubic angstroms








# Load flux data
data = np.loadtxt("heat_flux.dat")
timesteps = data[:, 0]
Jx = data[:, 1]
Jy = data[:, 2]
Jz = data[:, 3]



# TRUNCATE DATA IF NECESSARY (e.g if run for too long and kappa begins drifting)
# Define truncation time in picoseconds if slicing heat_flux
truncation_time = 1000 

# Convert truncation time to timesteps
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

# Multiply by V if already scaled by volume (code is for non-scaled fluxes)
# Jx = V * Jx
# Jy = V * Jy
# Jz = V * Jz



# Downsample flux data for plotting
timesteps_downsampled = timesteps[::s]
Jx_downsampled = Jx[::s]
Jy_downsampled = Jy[::s]
Jz_downsampled = Jz[::s]


# Plot raw and downsampled Jx
plt.figure(figsize=(10, 6))
plt.plot(timesteps, Jx, label="Raw Jx", alpha=0.5)
plt.plot(timesteps_downsampled, Jx_downsampled, label="Downsampled Jx", marker='o', linestyle='-', color='red', markersize=3, linewidth=0.5)
plt.xlabel("Timesteps")
plt.ylabel("Jx (Heat Flux)")
plt.title("Raw and Downsampled Jx")
plt.legend()
plt.grid(True)
plt.savefig("sampling.png")








# Open file for J0Jt.dat output
with open("J0Jt_python.dat", "w") as f_out:
    f_out.write("# Time-correlated data for Python Green-Kubo analysis\n")
    f_out.write("# Timestep Number-of-time-windows\n")
    f_out.write("# Index TimeDelta Ncount Jx*Jx Jy*Jy Jz*Jz\n")


# Compute ACF function with incremental updates
# More efficient version than the one implemented in green_kubo_heatmap.py
# This does not recalculate old windows already averaged
# green_kubo_heatmap.py version is fine since it only does it once
# That version is also more instructive/easier to understand


def compute_acf(J, p, s, d_timestep, C_vals, counts):

    # Compute the autocorrelation function C_J(Δt) at every d timestep.

    # J: The downsampled heat flux array.
    # p: Number of correlation samples.
    # s: Sample interval.
    # d_timestep: Current time step, must be a multiple of (p * s)
    # C_vals: numpy array of current C(0), C(s), C(2s), ..., C(d - s). len(C_vals) should equal p.
    # counts: numpy array of number of windows contributing to each corresponding C, is also of size p.


    for p_num in range(0, p):
        lag = p_num * s
        lag_downsampled = p_num

        # Find number of windows available for the current C(lag)
        n_datapoints = int(d_timestep / s)
        n_windows = n_datapoints - p_num - p + 1

        if n_windows > 0:
            old_counts = int(counts[p_num])
            new_counts = int(n_windows - old_counts)
            
            if new_counts > 0:
                C_vals_new_unscaled = 0.0
                
                for i in range(old_counts, n_windows):
                    loc0 = i
                    loct = loc0 + lag_downsampled
                    
                    if loct + p > len(J):
                        break  # Prevent index out of bounds error

                    J_0 = J[loc0:loc0 + p]
                    J_t = J[loct:loct + p]

                    C_vals_new_unscaled += (1/p) * np.dot(J_0, J_t)
       
                # Average new with old (i.e 'ave running' option in LAMMPS)
                C_vals[p_num] = ((C_vals[p_num] * old_counts) + C_vals_new_unscaled) / n_windows
                
                
                counts[p_num] = n_windows
    
    return C_vals, counts











# Scaling factor for Green-Kubo integral
scale = (convert / (kB * T * T * V)) * s * dt

# Storage for tracking conductivity over time
thermal_conductivity_x = []
thermal_conductivity_y = []
thermal_conductivity_z = []
time_evolution = []

total_timesteps = len(timesteps)

if d >= total_timesteps:
    print("ERROR: Not enough timesteps to process. Increase simulation length.")
else:
    C_vals_x = np.zeros(p)
    counts_x = np.zeros(p)
    
    C_vals_y = np.zeros(p)
    counts_y = np.zeros(p)
    
    C_vals_z = np.zeros(p)
    counts_z = np.zeros(p)

    for d_timestep in range(d, total_timesteps + 1, d):
        C_x, counts_x = compute_acf(Jx_downsampled, p, s, d_timestep, C_vals_x, counts_x)
        C_y, counts_y = compute_acf(Jy_downsampled, p, s, d_timestep, C_vals_y, counts_y)
        C_z, counts_z = compute_acf(Jz_downsampled, p, s, d_timestep, C_vals_z, counts_z)

        kappa_11 = np.trapz(C_x) * scale
        kappa_22 = np.trapz(C_y) * scale
        kappa_33 = np.trapz(C_z) * scale

        time_evolution.append(d_timestep * dt)
        thermal_conductivity_x.append(kappa_11)
        thermal_conductivity_y.append(kappa_22)
        thermal_conductivity_z.append(kappa_33)

        # Write to J0Jt.dat file
        with open("J0Jt_python.dat", "a") as f_out:
            f_out.write(f"{d_timestep - d} {p}\n")
            for k in range(p):
                time_delta = int(k * s)
                f_out.write(f"{k+1} {time_delta} {int(counts_x[k])} {C_x[k]:.6f} {C_y[k]:.6f} {C_z[k]:.6f}\n")

    # Compute average thermal conductivity
    if thermal_conductivity_x:
        kappa = (thermal_conductivity_x[-1] + thermal_conductivity_y[-1] + thermal_conductivity_z[-1]) / 3.0
        print(f"\nFinal Average Thermal Conductivity: {kappa:.2f} W/mK\n")
    else:
        print("ERROR: No valid conductivity values were computed.")



# Plot evolution of thermal conductivity only if data exists
if thermal_conductivity_x:
    plt.figure(figsize=(10, 6))
    plt.plot(time_evolution, thermal_conductivity_x, label=r"$k_{xx}$")
    plt.plot(time_evolution, thermal_conductivity_y, label=r"$k_{yy}$")
    plt.plot(time_evolution, thermal_conductivity_z, label=r"$k_{zz}$")
    plt.xlabel("Time (ps)")
    plt.ylabel("Thermal Conductivity (W/mK)")
    plt.title("Evolution of Thermal Conductivity Over Time")
    plt.legend()
    plt.grid(True)
    plt.savefig("kappa_evolution.png")
