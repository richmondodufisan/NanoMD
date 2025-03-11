import numpy as np
from ase.io import read, write
from ase.build import make_supercell
from ase import Atoms
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

# Load the unit cell from the CIF file
unit_cell = read('Si.cif')  # Replace with your CIF file path

# Desired box dimensions [x, y, z]
desired_box = np.array([100, 50, 50])

# Get the lattice vectors of the unit cell
cell = unit_cell.get_cell()

# Extract the lengths of the unit cell along x, y, z axes
unit_cell_lengths = np.linalg.norm(cell, axis=1)

# Calculate the number of times to replicate the unit cell along each axis
scaling_factors = np.ceil(desired_box / unit_cell_lengths).astype(int)
print("Scaling factors (number of repetitions along x, y, z):", scaling_factors)

# Create the supercell by replicating the unit cell according to the scaling factors
scaling_matrix = np.diag(scaling_factors)
supercell = make_supercell(unit_cell, scaling_matrix)






############## FOR COLORING VISUALIZATION  #############
# Define the coloring regions
x_min = 0.2 * desired_box[0]
x_max = 0.8 * desired_box[0]


lx = desired_box[0]
hot_left_min = 0
hot_left_max = 0+(lx*0.1)

cold_min = hot_left_max+(lx*0.35)
cold_max = cold_min+(lx*0.1)

hot_right_min = cold_max+(lx*0.35)
hot_right_max = hot_right_min+(lx*0.1)





# Lists to store atoms with their colors
trimmed_positions = []
trimmed_symbols = []
colors = []

for atom in supercell:
    pos = atom.position
    if (0 < pos[0] < desired_box[0]) and (0 < pos[1] < desired_box[1]) and (0 < pos[2] < desired_box[2]):
    
    
        trimmed_positions.append(pos)
        trimmed_symbols.append(atom.symbol)
        
        
        
        
        ############## FOR COLORING VISUALIZATION  #############
        
        # Assign colors based on position
        # HOT - COLD
        # if pos[0] < x_min:
            # colors.append('red')  # Leftmost region
        # elif pos[0] > x_max:
            # colors.append('blue')  # Rightmost region
        # else:
            # colors.append('gray')  # Default color
            
            
        # Assign colors based on position
        # HOT - COLD - HOT
        # if ((pos[0] > hot_left_min) and (pos[0] < hot_left_max))  or ((pos[0] > hot_right_min) and (pos[0] < hot_right_max)):
            # colors.append('red')  # Leftmost region
        # elif ((pos[0] > cold_min) and (pos[0] < cold_max)):
            # colors.append('blue')  # Rightmost region
        # else:
            # colors.append('gray')  # Default color
            
        
        # All one color
        colors.append('red')




# Create a new Atoms object with only the atoms inside the box
trimmed_supercell = Atoms(positions=trimmed_positions, symbols=trimmed_symbols, cell=desired_box)
trimmed_supercell.set_pbc(True)  # Enable periodic boundary conditions



############## FOR COLORING VISUALIZATION  #############
# Visualization using Matplotlib
fig, ax = plt.subplots(figsize=(10, 5))

# Increase radii and remove black boundaries
plot_atoms(trimmed_supercell, ax, colors=colors, radii=0.8)
ax.set_xticks([])
ax.set_yticks([])
ax.set_frame_on(False)  # Remove axis frame
plt.savefig("atom_config.png")
plt.show()




# Save the final trimmed structure
trimmed_supercell.write('Silicon_supercell_200x25x25.data', format='lammps-data')
trimmed_supercell.write('Silicon_supercell_200x25x25.xyz')
