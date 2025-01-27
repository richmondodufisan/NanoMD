import numpy as np
import matplotlib.pyplot as plt

def calculate_slope(filename):
    """
    Calculate the slope of temperature vs position for a given data file.
    This function will also handle variable chunks based on suffix.
    """
    # Load the data from the text file, skipping the first 4 lines
    data = np.loadtxt(filename, skiprows=1)

    # Extract the second and fourth columns
    x_position = data[:, 1]  # second column
    temperature = data[:, 3]  # fourth column

    # Determine the number of chunks and how many to skip based on suffix
    suffix = float(filename.split('_')[-1].replace('.txt', ''))
    n_chunks = int(20 * suffix)
    n_to_skip = int(0.2 * n_chunks) + 1

    # Slice data based on n_to_skip
    x_position = x_position[n_to_skip:-n_to_skip]
    temperature = temperature[n_to_skip:-n_to_skip]

    # Perform linear regression to fit a straight line (y = mx + b)
    slope, intercept = np.polyfit(x_position, temperature, 1)

    return slope

def plot_temperature_data(filename):
    """
    Plot the temperature data and best fit line for a single file.
    """
    # Load the data from the text file
    data = np.loadtxt(filename, skiprows=1)

    # Extract the second and fourth columns
    x_position = data[:, 1]
    temperature = data[:, 3]

    # Perform linear regression to fit a straight line
    suffix = float(filename.split('_')[-1].replace('.txt', ''))
    n_chunks = int(20 * suffix)
    n_to_skip = int(0.2 * n_chunks) + 1

    # Slice data based on n_to_skip
    x_position = x_position[n_to_skip:-n_to_skip]
    temperature = temperature[n_to_skip:-n_to_skip]

    slope, intercept = np.polyfit(x_position, temperature, 1)

    # Create the best-fit line
    best_fit_line = slope * x_position + intercept

    # Plot the data and the best-fit line
    plt.figure()
    plt.plot(x_position, temperature, 'o', label='Data points')  # Original data points
    plt.plot(x_position, best_fit_line, 'r-', label=f'Best fit line: y={slope:.2f}x+{intercept:.2f}')  # Best-fit line
    plt.xlabel('X Position')
    plt.ylabel('Temperature')
    plt.title('X Position vs Temperature with Best Fit Line')
    plt.legend()
    plt.grid(True)

    plt.savefig(f'{filename.split("/")[-1].replace(".txt", "")}_temp_profile.png')

    # Show the plot
    plt.show()

    # Print the slope of the line
    print(f"The slope of the best-fit line is: {slope:.2f}")

# Main code block for single-file plotting (for testing)
if __name__ == "__main__":
    plot_temperature_data('./lammps_out/temp_profile_1.5.txt')
