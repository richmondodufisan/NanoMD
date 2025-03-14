import matplotlib.pyplot as plt






# File name (change this to your actual file)
file_name = "equilibrium.log"

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
plt.xlabel("Time (femtoseconds)")
plt.ylabel("Temperature (K)")
plt.title("Timestep vs Temperature - Equilibration Phase")
plt.grid(True)
plt.savefig("equilibration.png")
plt.show()






# # File name (change this to your actual file)
# file_name = "thermal_conductivity.log"

# # Read data from file
# timesteps = []
# temperatures = []

# kappa_11 = []
# kappa_22 = []
# kappa_33 = []

# with open(file_name, "r") as file:
    # for line in file:
        # if line.strip():  # Ignore empty lines
            # parts = line.split()  # Split by whitespace (change to .split(',') if CSV)
            # try:
                # timesteps.append(float(parts[0]))
                # temperatures.append(float(parts[1]))
                
                # kappa_11.append(float(parts[5]))
                # kappa_22.append(float(parts[6]))
                # kappa_33.append(float(parts[7]))
            # except ValueError:
                # print(f"Skipping invalid line: {line.strip()}")

# # Plot the data
# plt.figure(figsize=(8, 5))
# plt.plot(timesteps, temperatures)
# plt.xlabel("Time (femtoseconds)")
# plt.ylabel("Temperature (K)")
# plt.title("Timestep vs Temperature - Constant NVE Phase")
# plt.grid(True)
# plt.savefig("NVE.png")
# plt.show()


# # Plot the data
# plt.figure(figsize=(8, 5))
# plt.plot(timesteps, kappa_11, label="κ_11")
# plt.plot(timesteps, kappa_22, label="κ_22")
# plt.plot(timesteps, kappa_33, label="κ_33")
# plt.xlabel("Time (femtoseconds)")
# plt.ylabel("Thermal Conductivity (W/(m·K))")
# plt.title("Timestep vs Thermal Conductivity")
# plt.grid(True)
# plt.legend()  # Add legend for kappa values
# plt.savefig("time_v_kappa.png")
# plt.show()