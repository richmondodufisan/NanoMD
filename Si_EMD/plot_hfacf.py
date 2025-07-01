import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd

# === STEP 1: Read the file ===
with open("J0Jt_python.dat", "r") as f:
    lines = f.readlines()

# === STEP 2: Parse into blocks ===
blocks = []
current_block = {"header": None, "data": []}

for line in lines[3:]:  # Skip the first 3 header lines
    if re.match(r"^\d+\s+\d+\s*$", line):
        if current_block["header"] and current_block["data"]:
            blocks.append(current_block)
        timestep, num_windows = map(int, line.strip().split())
        current_block = {"header": timestep, "data": []}
    else:
        current_block["data"].append(list(map(float, line.strip().split())))

if current_block["header"] and current_block["data"]:
    blocks.append(current_block)

# === STEP 3: Compute integrals per block ===
real_times = []
integrals = []

for block in blocks:
    timestep = block["header"]
    data = np.array(block["data"])
    if data.shape[1] < 6:
        continue  # Skip malformed rows

    # Extract Cx, Cy, Cz
    Cx = data[:, 3]
    Cy = data[:, 4]
    Cz = data[:, 5]
    avg_C = (Cx + Cy + Cz) / 3.0

    # Integrate avg_C over time deltas
    time_deltas = data[:, 1]
    integral = np.trapz(avg_C, time_deltas)

    real_times.append(timestep)
    integrals.append(integral)

# === STEP 4: Normalize ===
# Normalize by first value
# normalized_integrals = [val / integrals[0] for val in integrals]

# OR: Normalize by maximum value
normalized_integrals = [val / max(integrals) for val in integrals]

# === STEP 5: Plot results (Journal Quality) ===
plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "lines.linewidth": 2,
    "lines.markersize": 6
})

plt.figure(figsize=(10, 6))
plt.plot(real_times, normalized_integrals, marker='o', linestyle='-')
plt.title("Normalized Integral of Averaged HFACF per Interval")
plt.xlabel("Real Time (ps)")
plt.ylabel("Normalized Integrated ACF")
plt.grid(True)
plt.tight_layout()
plt.savefig("normalized_integrated_acf_vs_time.png", dpi=600)
# plt.show()


# Optional: save results to CSV
df = pd.DataFrame({
    "Timestep Origin": real_times,
    "Normalized Integrated ACF": normalized_integrals
})
df.to_csv("normalized_acf_integrals.csv", index=False)
