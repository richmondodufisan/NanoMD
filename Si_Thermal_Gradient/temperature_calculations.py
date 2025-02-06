import numpy as np
import matplotlib.pyplot as plt

def calculate_slope(filename):
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
    n_chunks = int(40 * suffix)
    
    n_to_skip_front = int(0.075 * n_chunks)
    n_to_skip_back = int(0.5 * n_chunks) + n_to_skip_front

    # Slice data based on n_to_skip
    x_position_fit = x_position[n_to_skip_front:-n_to_skip_back]
    temperature_fit = temperature[n_to_skip_front:-n_to_skip_back]
    
    # x_position_fit = x_position[7:-43]
    # temperature_fit = temperature[7:-43]

    slope, intercept = np.polyfit(x_position_fit, temperature_fit, 1)

    # Create the best-fit line
    best_fit_line = slope * x_position_fit + intercept

    # Plot the data and the best-fit line
    plt.figure()
    plt.plot(x_position, temperature, 'o', label='Data points')  # Original data points
    plt.plot(x_position_fit, best_fit_line, 'r-', label=f'Best fit line: y={slope:.2f}x+{intercept:.2f}')  # Best-fit line
    plt.xlabel('X Position')
    plt.ylabel('Temperature')
    plt.title('X Position vs Temperature with Best Fit Line')
    plt.legend()
    plt.grid(True)

    plt.savefig(f'{filename.split("/")[-1].replace(".txt", "")}_temp_profile_2.0.png')

    # Show the plot
    plt.show()

    # Print the slope of the line
    print(f"The slope of the best-fit line is: {slope:.2f}")
    
    return slope

# Main code block for single-file plotting (for testing)
if __name__ == "__main__":
    plot_temperature_data('./output_new/temp_profile_2.0.txt')
