from heat_flux_calculations import read_heat_flux_data, savitzky_golay_smooth, get_steady_state_flux
import temperature_calculations

import matplotlib.pyplot as plt
import numpy as np

import pdb

# Function to calculate box volume based on the suffix
# The middle box of the "1.0" simulation is 120 x 50 x 50
def calculate_box_volume(suffix):
    """
    Calculate the box volume by scaling the base dimensions (120 x 50 x 50)
    with the suffix.
    """
    base_length = 120
    base_width = 50
    base_height = 50
    return (base_length * suffix) * (base_width * suffix) * (base_height * suffix)


extensive_flux = [776.5761663468576, 940.1539851246981, 1430.1208944273103]
steady_state_fluxes = []

suffixes = [0.6, 0.7, 0.8]
box_volumes = []


for suffix in suffixes:
  box_volumes.append(calculate_box_volume(suffix))
  
for flux, volume in zip(extensive_flux, box_volumes):
  steady_state_fluxes.append(flux/volume)

# Plot scaled steady-state flux vs box volume
plt.figure(figsize=(10, 6))
plt.plot(box_volumes, steady_state_fluxes, marker='o', linestyle='-', color='blue', label='Scaled Steady-State Flux')

plt.xlabel('Box Volume')
plt.ylabel('Scaled Steady-State Heat Flux')
plt.title('Scaled Steady-State Heat Flux vs Box Volume')
plt.grid(True)
plt.legend()
plt.savefig('heat_flux_vs_volume')
plt.show()


slopes = [-10.20, -9.06, -8.0]

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



# Calculate and plot thermal conductivity vs box volume
plot_thermal_conductivity_vs_volume(steady_state_fluxes, slopes, box_volumes)
