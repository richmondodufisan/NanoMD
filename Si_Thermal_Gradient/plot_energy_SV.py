import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

# Function to read and parse the data from the text file
def read_energy_data(filename):
    time = []
    kinetic_energy = []
    potential_energy = []
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            # Split the line by comma and extract kinetic and potential energy
            kinetic, potential = map(float, line.strip().split(','))
            # The time is based on 100 nanoseconds increments
            time_value = i * 100  # each entry represents 100 nanoseconds
            time.append(time_value)
            kinetic_energy.append(kinetic)
            potential_energy.append(potential)
    return time, kinetic_energy, potential_energy

# Function to apply Savitzky-Golay smoothing
def savitzky_golay_smooth(data, window_size, poly_order):
    smoothed_data = savgol_filter(data, window_size, poly_order)
    return smoothed_data

# Function to plot the original and smoothed data
def plot_data(time, original_data, smoothed_data, data_label):
    plt.figure(figsize=(10, 6))
    
    # Plot original data
    plt.plot(time, original_data, label=f"Original {data_label}", alpha=0.5)
    
    # Plot smoothed data using Savitzky-Golay filter
    plt.plot(time, smoothed_data, label=f"Savitzky-Golay Smoothed {data_label}", zorder=2)
    
    plt.xlabel('Time (nanoseconds)')
    plt.ylabel(data_label)
    plt.title(f'{data_label} over Time (Savitzky-Golay Smoothing)')
    plt.grid(True)
    plt.legend()
    plt.show()

# Main code
filename = './Data/200X50X50/seed_1/energy_output_middle.txt'  # Replace with the path to your file
time, kinetic_energy, potential_energy = read_energy_data(filename)

# Set parameters for Savitzky-Golay filter
window_size = 501  # Needs to be odd and adapted based on your data length
poly_order = 3

# Apply Savitzky-Golay smoothing to kinetic and potential energy
smoothed_ke_savgol = savitzky_golay_smooth(kinetic_energy, window_size, poly_order)
smoothed_pe_savgol = savitzky_golay_smooth(potential_energy, window_size, poly_order)

# Calculate total energy (kinetic + potential)
total_energy = np.array(kinetic_energy) + np.array(potential_energy)
smoothed_total_energy = smoothed_ke_savgol + smoothed_pe_savgol

# Plot the original and smoothed kinetic energy data
plot_data(time, kinetic_energy, smoothed_ke_savgol, "Kinetic Energy")

# Plot the original and smoothed potential energy data
plot_data(time, potential_energy, smoothed_pe_savgol, "Potential Energy")

# Plot the total energy (kinetic + potential) and its smoothed version
plot_data(time, total_energy, smoothed_total_energy, "Total Energy")

# Print the steady-state kinetic, potential, and total energy as the last values from the smoothed data
steady_state_ke_savgol = smoothed_ke_savgol[-1]
steady_state_pe_savgol = smoothed_pe_savgol[-1]
steady_state_total_savgol = smoothed_total_energy[-1]

print(f"Steady state kinetic energy (Savitzky-Golay): {steady_state_ke_savgol}")
print(f"Steady state potential energy (Savitzky-Golay): {steady_state_pe_savgol}")
print(f"Steady state total energy (Savitzky-Golay): {steady_state_total_savgol}")
