import numpy as np
from ase.io import read, write
from ase.build import make_supercell
from ase import Atoms
from ase.visualize import view
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt



# Load the unit cell from the CIF file
unit_cell = read('Si.cif')  # Replace with your CIF file path

# Desired final box dimensions [x, y, z]
desired_box = np.array([200, 50, 50])

# Rotation angles (in degrees)
left_rotation_angle = 0 
right_rotation_angle = 15 

# Split the desired box in half along the x-axis
left_box_dimensions = desired_box / np.array([2, 1, 1])
right_box_dimensions = desired_box / np.array([2, 1, 1])

# Get the lattice vectors of the unit cell
cell = unit_cell.get_cell()
unit_cell_lengths = np.linalg.norm(cell, axis=1)




# Function to create a rotated and replicated supercell
def create_rotated_supercell(unit_cell, box_dimensions, rotation_angle):
    # Rotate the unit cell
    rotated_cell = unit_cell.copy()
    rotated_cell.rotate(rotation_angle, 'z', center='COM')
    rotated_cell.center()  # Ensure it's centered after rotation

    # Add buffer to ensure complete trimming after rotation
    buffer = np.array([5.0, 5.0, 5.0])
    scaling_factors = np.ceil((box_dimensions + buffer) / unit_cell_lengths).astype(int)
    scaling_matrix = np.diag(scaling_factors)
    supercell = make_supercell(rotated_cell, scaling_matrix)

    # Trim atoms outside the desired box
    trimmed_positions = []
    trimmed_symbols = []
    for atom in supercell:
        pos = atom.position
        if (0 <= pos[0] < box_dimensions[0]) and (0 <= pos[1] < box_dimensions[1]) and (0 <= pos[2] < box_dimensions[2]):
            trimmed_positions.append(pos)
            trimmed_symbols.append(atom.symbol)

    trimmed_supercell = Atoms(positions=trimmed_positions, symbols=trimmed_symbols, cell=box_dimensions)
    trimmed_supercell.set_pbc(True)
    return trimmed_supercell



# Create the left and right supercells
left_supercell = create_rotated_supercell(unit_cell, left_box_dimensions, left_rotation_angle)
right_supercell = create_rotated_supercell(unit_cell, right_box_dimensions, -right_rotation_angle)



# Print atom counts for verification
print(f"Left atoms:  {len(left_supercell)}")
print(f"Right atoms: {len(right_supercell)}")



# Shift the right supercell to join with the left one
right_supercell.positions[:, 0] += left_box_dimensions[0]

# Combine the two supercells
combined_supercell = left_supercell + right_supercell



# Final trimming to ensure box boundaries are respected
final_positions = []
final_symbols = []
colors = []


for atom in combined_supercell:
    
    pos = atom.position
    if (0 <= pos[0] < desired_box[0]) and (0 <= pos[1] < desired_box[1]) and (0 <= pos[2] < desired_box[2]):
        final_positions.append(pos)
        final_symbols.append(atom.symbol)
        
        
        # Define GB center and width (along x-axis)
        gb_center = desired_box[0] / 2
        gb_half_width = 5  # adjust this to change thickness
        
        right = gb_center + gb_half_width
        left = gb_center - gb_half_width
        
        ############## FOR COLORING VISUALIZATION  #############
        
        # # Assign colors based on position
        if ((pos[0] > left) and (pos[0] < right)):
            colors.append((1.0, 1.0, 1.0))  # white
        else:
            colors.append((0.4, 0.4, 0.4))  # dark gray (default color)
        
        
        
final_supercell = Atoms(positions=final_positions, symbols=final_symbols, cell=desired_box)
final_supercell.set_pbc(True)


############# FOR COLORING VISUALIZATION  #############
# Visualization using Matplotlib
fig, ax = plt.subplots(figsize=(10, 5))

# Increase radii and remove black boundaries
plot_atoms(final_supercell, ax, colors=colors, radii=1.5)
ax.set_xticks([])
ax.set_yticks([])
ax.set_frame_on(False)  # Remove axis frame
plt.savefig("atom_config.png")
plt.show()



# Visualize and save
# view(final_supercell)
final_supercell.write('Silicon_grain_boundary_50x15x15.data', format='lammps-data')
final_supercell.write('Silicon_grain_boundary_50x15x15.xyz')
