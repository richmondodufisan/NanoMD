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
data = np.loadtxt("test_flux.dat")
timesteps = data[:, 0]
Jx = data[:, 1]
Jy = data[:, 2]
Jz = data[:, 3]

# Multiply by V if already scaled by volume (code is for non-scaled fluxes)
Jx = V * Jx
Jy = V * Jy
Jz = V * Jz

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
plt.savefig("sampling_old.png")
# plt.show()

# Open file for J0Jt.dat output
with open("J0Jt_python_old.dat", "w") as f_out:

# Index TimeDelta Ncount c_flux[1]*c_flux[1] c_flux[2]*c_flux[2] c_flux[3]*c_flux[3]
    f_out.write("# Time-correlated data for Python Green-Kubo analysis\n")
    f_out.write("# Timestep Number-of-time-windows\n")
    f_out.write("# Index TimeDelta Ncount Jx*Jx Jy*Jy Jz*Jz\n")



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

# Scaling factor for Green-Kubo integral
# if flux is NOT scaled by volume
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
    for d_timestep in range(d, total_timesteps + 1, d):
        C_x, counts_x = compute_acf(Jx_downsampled, p, s, d_timestep)
        C_y, counts_y = compute_acf(Jy_downsampled, p, s, d_timestep)
        C_z, counts_z = compute_acf(Jz_downsampled, p, s, d_timestep)

        kappa_11 = np.trapz(C_x) * scale
        kappa_22 = np.trapz(C_y) * scale
        kappa_33 = np.trapz(C_z) * scale

        time_evolution.append(d_timestep * dt)
        thermal_conductivity_x.append(kappa_11)
        thermal_conductivity_y.append(kappa_22)
        thermal_conductivity_z.append(kappa_33)

        # Write to J0Jt.dat file
        with open("J0Jt_python_old.dat", "a") as f_out:
            f_out.write(f"{d_timestep - d} {p}\n")
            
            for k in range(p):
                time_delta = int(k * s)
                f_out.write(f"{k+1} {time_delta} {int(counts_x[k])} {C_x[k]:.6f} {C_y[k]:.6f} {C_z[k]:.6f}\n")

    # Compute average thermal conductivity only if data exists
    if thermal_conductivity_x:
        kappa = (thermal_conductivity_x[-1] + thermal_conductivity_y[-1] + thermal_conductivity_z[-1]) / 3.0
        print(f"\nFinal Average Thermal Conductivity: {kappa:.2f} W/mK\n")
    else:
        print("ERROR: No valid conductivity values were computed.")

# Plot evolution of thermal conductivity
if thermal_conductivity_x:
    plt.figure(figsize=(10, 6))
    plt.plot(time_evolution, thermal_conductivity_x, label=r"$k_{xx}$", linestyle='-', marker='o')
    plt.plot(time_evolution, thermal_conductivity_y, label=r"$k_{yy}$", linestyle='-', marker='s')
    plt.plot(time_evolution, thermal_conductivity_z, label=r"$k_{zz}$", linestyle='-', marker='^')
    plt.xlabel("Time (ps)")
    plt.ylabel("Thermal Conductivity (W/mK)")
    plt.title("Evolution of Thermal Conductivity Over Time")
    plt.legend()
    plt.grid(True)
    plt.savefig("kappa_evolution_old.png")
    # plt.show()
