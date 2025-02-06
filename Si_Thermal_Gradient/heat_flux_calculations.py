import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

# Function to read and parse the data from the text file
def read_heat_flux_data(filename):
    time = []
    heat_flux = []
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if i == 0:  # Skip the first line
                continue
            # Assume each subsequent line contains only the heat flux value
            flux_value = float(line.strip())
            # The time is based on 100 nanoseconds increments
            time_value = (i - 1) * 100  # Subtract 1 to account for the skipped first line
            time.append(time_value)
            heat_flux.append(flux_value)
    return time, heat_flux


# Function to apply Savitzky-Golay smoothing
def savitzky_golay_smooth(time, heat_flux, window_size, poly_order):
    smoothed_flux = savgol_filter(heat_flux, window_size, poly_order)
    return smoothed_flux


# Function to get steady state heat flux (last value of smoothed data)
def get_steady_state_flux(smoothed_flux):
    return smoothed_flux[-1]


# Function to plot the heat flux data along with the smoothed data
def plot_heat_flux(time, heat_flux, smoothed_flux, output_file=None):
    plt.figure(figsize=(10, 6))
    
    # Plot original data
    plt.plot(time, heat_flux, label="Original Heat Flux", alpha=0.5)
    
    # Plot smoothed data using Savitzky-Golay filter
    plt.plot(time, smoothed_flux, color='green', label="Savitzky-Golay Smoothed Heat Flux", zorder=2)
    
    plt.xlabel('Time (nanoseconds)')
    plt.ylabel('Heat Flux')
    plt.title('Heat Flux in x-direction over Time (Savitzky-Golay Smoothing)')
    plt.grid(True)
    plt.legend()
    
    if output_file:
        plt.savefig(output_file)
    plt.show()


# Main code block for single-file processing
if __name__ == "__main__":
    filename = './lammps_current/heat_flux_output_middle_left_2.0.txt'  # Replace with the path to your file
    time, heat_flux = read_heat_flux_data(filename)

    # Set parameters for Savitzky-Golay filter
    window_size = 5001  # Needs to be odd
    poly_order = 3

    # Apply Savitzky-Golay smoothing
    smoothed_flux = savitzky_golay_smooth(time, heat_flux, window_size, poly_order)


    # Plot the original and smoothed data
    plot_heat_flux(time, heat_flux, smoothed_flux, 'flux_v_time2.0.png')

    # Print the steady-state heat flux
    steady_state_flux = get_steady_state_flux(smoothed_flux)
    print(f"Steady state heat flux (Savitzky-Golay): {steady_state_flux}")
