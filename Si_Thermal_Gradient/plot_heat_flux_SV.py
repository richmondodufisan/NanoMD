import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

# Function to read and parse the data from the text file
def read_heat_flux_data(filename):
    time = []
    heat_flux = []
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if "Heat flux in x-direction for middle region:" in line:
                # Extract the heat flux value
                flux_value = float(line.split(":")[1].strip())
                # The time is based on 100 nanoseconds increments
                time_value = i * 100  # each entry represents 100 nanoseconds
                time.append(time_value)
                heat_flux.append(flux_value)
    return time, heat_flux

# Function to apply Savitzky-Golay smoothing
def savitzky_golay_smooth(time, heat_flux, window_size, poly_order):
    smoothed_flux = savgol_filter(heat_flux, window_size, poly_order)
    return smoothed_flux

# Function to plot the heat flux data along with the smoothed data
def plot_heat_flux(time, heat_flux, smoothed_flux_savgol):
    plt.figure(figsize=(10, 6))
    
    # Plot original data
    plt.plot(time, heat_flux, label="Original Heat Flux", alpha=0.5)
    
    # Plot smoothed data using Savitzky-Golay filter
    plt.plot(time, smoothed_flux_savgol, color='green', label="Savitzky-Golay Smoothed Heat Flux", zorder=2)
    
    plt.xlabel('Time (nanoseconds)')
    plt.ylabel('Heat Flux')
    plt.title('Heat Flux in x-direction over Time (Savitzky-Golay Smoothing)')
    plt.grid(True)
    plt.legend()
    plt.show()

# Main code
filename = './Data/200X50X50/seed_1/heat_flux_output_middle.txt'  # Replace with the path to your file
time, heat_flux = read_heat_flux_data(filename)

# Set parameters for Savitzky-Golay filter
window_size = 5001  # Needs to be odd
poly_order = 3

# Apply Savitzky-Golay smoothing
smoothed_flux_savgol = savitzky_golay_smooth(time, heat_flux, window_size, poly_order)

# Plot the original and smoothed data
plot_heat_flux(time, heat_flux, smoothed_flux_savgol)

# Print the steady state heat flux as the last value from the smoothed data
steady_state_heat_flux_savgol = smoothed_flux_savgol[-1]
print(f"Steady state heat flux (Savitzky-Golay): {steady_state_heat_flux_savgol}")
