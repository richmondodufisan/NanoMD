import numpy as np
import matplotlib.pyplot as plt

# Load the data from the text file, skipping the first 4 lines
data = np.loadtxt('Data/200X50X50/temp.txt', skiprows=4)

# Extract the second and fourth columns
x_position = data[:, 1]  # second column
temperature = data[:, 3]  # fourth column

# Perform linear regression to fit a straight line (y = mx + b)
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

# Show the plot
plt.show()

# Print the slope of the line
print(f"The slope of the best-fit line is: {slope:.2f}")
