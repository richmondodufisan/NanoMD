import matplotlib.pyplot as plt
import numpy as np

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

# Function to smooth the heat flux data by chunking into bins and averaging
def smooth_heat_flux(time, heat_flux, bin_size_ns):
    smoothed_time = []
    smoothed_flux = []
    
    start_time = time[0]
    current_time_bin_start = start_time
    flux_sum = 0
    count = 0
    
    # Iterate through the data and chunk by bin size (based on time)
    for t, flux in zip(time, heat_flux):
        if t < current_time_bin_start + bin_size_ns:
            flux_sum += flux
            count += 1
        else:
            # Bin is full, calculate average
            avg_flux = flux_sum / count if count > 0 else 0
            avg_time = current_time_bin_start + bin_size_ns / 2  # Center of the bin
            smoothed_time.append(avg_time)
            smoothed_flux.append(avg_flux)
            
            # Reset for the next bin
            current_time_bin_start += bin_size_ns
            flux_sum = flux
            count = 1  # reset count with current data point

    # Handle the last bin if there are remaining data points
    if count > 0:
        avg_flux = flux_sum / count
        avg_time = current_time_bin_start + bin_size_ns / 2
        smoothed_time.append(avg_time)
        smoothed_flux.append(avg_flux)
    
    return smoothed_time, smoothed_flux

# Function to plot the heat flux data along with the smoothed data
def plot_heat_flux(time, heat_flux, smoothed_time, smoothed_flux):
    plt.figure(figsize=(10, 6))
    
    # Plot original data
    plt.plot(time, heat_flux, label="Original Heat Flux", alpha=0.5)
    
    # Plot smoothed data as dots and connect them with lines
    plt.scatter(smoothed_time, smoothed_flux, color='red', label="Smoothed Heat Flux", zorder=3)
    plt.plot(smoothed_time, smoothed_flux, color='red', zorder=2)
    
    plt.xlabel('Time (nanoseconds)')
    plt.ylabel('Heat Flux')
    plt.title('Heat Flux in x-direction over Time')
    plt.grid(True)
    plt.legend()
    plt.show()

# Main code
filename = './Data/400X100X100/seed_1/heat_flux_output_middle.txt'  # Replace with the path to your file
time, heat_flux = read_heat_flux_data(filename)

# time = time[:-35000]
# heat_flux = heat_flux[:-35000]

# Set the bin size in nanoseconds (e.g., 10,000 nanoseconds)
bin_size_ns = 10000

# Smooth the data based on the bin size
smoothed_time, smoothed_flux = smooth_heat_flux(time, heat_flux, bin_size_ns)

# Plot the original and smoothed data
plot_heat_flux(time, heat_flux, smoothed_time, smoothed_flux)

# Print the steady state heat flux as the last value from the smoothed data
steady_state_heat_flux = smoothed_flux[-1]
print(f"Steady state heat flux: {steady_state_heat_flux}")
