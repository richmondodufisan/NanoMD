import matplotlib.pyplot as plt

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

# Function to plot the heat flux data
def plot_heat_flux(time, heat_flux):
    plt.figure(figsize=(10, 6))
    plt.plot(time, heat_flux, label="Heat Flux in x-direction")
    plt.xlabel('Time (nanoseconds)')
    plt.ylabel('Heat Flux')
    plt.title('Heat Flux in x-direction over Time')
    plt.grid(True)
    plt.legend()
    plt.show()

# Main code
filename = 'heat_flux_output_middle.txt'  # Replace with the path to your file
time, heat_flux = read_heat_flux_data(filename)
plot_heat_flux(time, heat_flux)
