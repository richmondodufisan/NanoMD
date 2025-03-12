import matplotlib.pyplot as plt

# File name (change this to your actual file)
file_name = "data.txt"

# Read data from file
timesteps = []
temperatures = []

with open(file_name, "r") as file:
    for line in file:
        if line.strip():  # Ignore empty lines
            parts = line.split()  # Split by whitespace (change to .split(',') if CSV)
            try:
                timesteps.append(float(parts[0]))
                temperatures.append(float(parts[1]))
            except ValueError:
                print(f"Skipping invalid line: {line.strip()}")

# Plot the data
plt.figure(figsize=(8, 5))
plt.plot(timesteps, temperatures, marker='o', linestyle='-')
plt.xlabel("Timestep")
plt.ylabel("Temperature")
plt.title("Timestep vs Temperature")
plt.grid(True)
plt.show()
