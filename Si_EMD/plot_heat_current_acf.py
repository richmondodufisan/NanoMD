import numpy as np
import matplotlib.pyplot as plt
import re

# === STEP 1: Read the file ===
with open("J0Jt_python.dat", "r") as f:
    lines = f.readlines()

# === STEP 2: Skip headers and organize data into blocks ===
# Each block starts with a line like "0 100" or "2000 100"
blocks = []
current_block = {"header": None, "data": []}

for line in lines[3:]:  # Skip the first 3 lines (text headers)
    if re.match(r"^\d+\s+\d+\s*$", line):
        # If we hit a new header line, save the current block (if it has data)
        if current_block["header"] and current_block["data"]:
            blocks.append(current_block)
        # Start a new block with this header (time_origin, num_windows)
        timestep, num_windows = map(int, line.strip().split())
        current_block = {"header": (timestep, num_windows), "data": []}
    else:
        # Add correlation data row to the current block
        current_block["data"].append(line.strip().split())

# Don't forget to append the last block!
if current_block["header"] and current_block["data"]:
    blocks.append(current_block)

# === STEP 3: Accumulate and average ACF values by time shift (TimeDelta) ===
acf_accumulator = {}  # Dictionary to sum ACFs per shift
acf_counts = {}       # Dictionary to count how many values are summed

# Loop through each block (each time origin)
for block in blocks:
    for row in block["data"]:
        _, timedelta, ncount, jx, jy, jz = row
        timedelta = int(timedelta)
        ncount = int(ncount)
        
        # Skip entries with 0 contributing samples
        if ncount == 0:
            continue
        
        # Compute the average of the x, y, z flux correlations
        avg_acf = (float(jx) + float(jy) + float(jz)) / 3.0
        
        # Accumulate total value and count for this time shift
        if timedelta not in acf_accumulator:
            acf_accumulator[timedelta] = 0.0
            acf_counts[timedelta] = 0
        
        acf_accumulator[timedelta] += avg_acf
        acf_counts[timedelta] += 1

# === STEP 4: Compute final averaged ACF for each time shift ===
timedeltas = sorted(acf_accumulator.keys())
acf_avg = [acf_accumulator[t] / acf_counts[t] for t in timedeltas]

# === STEP 5: Plot the result ===
plt.figure(figsize=(10, 6))
plt.plot(timedeltas, acf_avg, marker='o', linestyle='-')
plt.title("Averaged Heat Flux/Current Autocorrelation Function (HCACF)")
plt.xlabel("Time Shift (Δt)")  # This is the lag/shift time, not simulation time
plt.ylabel("Average ACF")
plt.grid(True)
plt.tight_layout()
plt.savefig('hcacf_plot.png')
plt.show()
