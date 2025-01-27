from heat_flux_calculations import read_heat_flux_data, savitzky_golay_smooth, get_steady_state_flux
import temperature_calculations

import matplotlib.pyplot as plt
import numpy as np

import pdb

# Function to calculate box volume based on the suffix
def calculate_box_volume(suffix):
    """
    Calculate the box volume by scaling the base dimensions (120 x 50 x 50)
    with the suffix.
    """
    base_length = 120
    base_width = 50
    base_height = 50
    return (base_length * suffix) * (base_width * suffix) * (base_height * suffix)


# Function to process files and plot steady-state flux vs. box volume
def plot_steady_state_vs_volume(filenames, output_plot='steady_state_vs_volume.png'):
    steady_state_fluxes = []
    box_volumes = []

    for filename in filenames:
        # Extract suffix from filename (e.g., "1.0" from "heat_flux_output_middle_1.0.txt")
        suffix = float(filename.split('_')[-1].replace('.txt', ''))

        # Calculate the box volume
        box_volume = calculate_box_volume(suffix)

        # Read and process data
        time, heat_flux = read_heat_flux_data(filename)
        smoothed_flux = savitzky_golay_smooth(time, heat_flux, WINDOW_SIZE, POLY_ORDER)
        steady_state_flux = get_steady_state_flux(smoothed_flux)

        # Scale steady-state heat flux by box volume
        scaled_flux = steady_state_flux / box_volume

        # Append results for plotting
        box_volumes.append(box_volume)
        steady_state_fluxes.append(scaled_flux)

    # Plot scaled steady-state flux vs box volume
    plt.figure(figsize=(10, 6))
    plt.plot(box_volumes, steady_state_fluxes, marker='o', linestyle='-', color='blue', label='Scaled Steady-State Flux')

    plt.xlabel('Box Volume')
    plt.ylabel('Scaled Steady-State Heat Flux')
    plt.title('Scaled Steady-State Heat Flux vs Box Volume')
    plt.grid(True)
    plt.legend()
    plt.savefig(output_plot)
    plt.show()

    return box_volumes, steady_state_fluxes


# Function to plot temperature slope vs. box volume
def plot_temperature_slope_vs_volume():
    slopes = []
    box_volumes = []

    # Exact suffixes to avoid floating-point approximation errors
    suffixes = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

    for suffix in suffixes:
        # Generate the filename for the temperature profile
        filename = f"./lammps_out/temp_profile_{suffix}.txt"

        try:
            # Calculate the temperature slope using the helper file
            slope = temperature_calculations.calculate_slope(filename)

            # Calculate the box volume
            box_volume = calculate_box_volume(suffix)

            # Append results for plotting
            slopes.append(slope)
            box_volumes.append(box_volume)

        except FileNotFoundError:
            print(f"File {filename} not found, skipping.")

    # Plot temperature slope vs box volume
    plt.figure(figsize=(10, 6))
    plt.plot(box_volumes, slopes, marker='o', linestyle='-', color='red', label='Temperature Slope')

    plt.xlabel('Box Volume')
    plt.ylabel('Temperature Slope')
    plt.title('Temperature Slope vs Box Volume')
    plt.grid(True)
    plt.legend()
    plt.savefig('temperature_slope_vs_volume.png')
    plt.show()

    return slopes, box_volumes


# Function to calculate thermal conductivity and plot it
def plot_thermal_conductivity_vs_volume(steady_state_fluxes, slopes, box_volumes):
    thermal_conductivities = []

    # Calculate thermal conductivity for each data point
    for flux, slope, volume in zip(steady_state_fluxes, slopes, box_volumes):
    
        # Calculate temperature gradient (temperature slope / box size)
        temperature_gradient = slope

        # Calculate thermal conductivity with the negative sign
        thermal_conductivity = -flux / temperature_gradient
        
        # Convert to W/(m.K)
        thermal_conductivity = thermal_conductivity * ((1.60218e-19)/(1e-22))
        # pdb.set_trace()
        
        thermal_conductivities.append(thermal_conductivity)

    # Plot thermal conductivity vs box volume
    plt.figure(figsize=(10, 6))
    plt.plot(box_volumes, thermal_conductivities, marker='o', linestyle='-', color='green', label='Thermal Conductivity')

    plt.xlabel('Box Volume (cubic angstroms)')
    plt.ylabel('Thermal Conductivity (W/(m.K)')
    plt.title('Thermal Conductivity vs Box Volume')
    plt.grid(True)
    plt.legend()
    plt.savefig('thermal_conductivity_vs_volume.png')
    plt.show()



# Main code block
if __name__ == "__main__":
    # Parameters for Savitzky-Golay smoothing
    WINDOW_SIZE = 5001  # Needs to be odd
    POLY_ORDER = 3

    # Generate filenames based on exact suffixes
    suffixes = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    filenames = [f"./lammps_out/heat_flux_output_middle_{suffix}.txt" for suffix in suffixes]

    # Process files and plot the steady-state flux vs box volume
    box_volumes, steady_state_fluxes = plot_steady_state_vs_volume(filenames)

    # Plot the temperature slope vs box volume
    slopes, _ = plot_temperature_slope_vs_volume()

    # Calculate and plot thermal conductivity vs box volume
    plot_thermal_conductivity_vs_volume(steady_state_fluxes, slopes, box_volumes)
