import numpy as np
from ase.io import read
from ase.build import make_supercell
from ase import Atoms
from ase.visualize import view

# Load the unit cell from the CIF file
unit_cell = read('Si.cif')  # Replace with your CIF file path

# Desired box dimensions [x, y, z]
desired_box = np.array([50, 25, 25])

# Rotation angles (in degrees)
left_rotation_angle = 15  # Clockwise rotation for the left side
right_rotation_angle = 15  # Clockwise rotation for the right side (negative for counter-clockwise in this setup)

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
    rotated_cell.rotate(rotation_angle, 'z', center='COM') # Rotate about center of mass for better results.

    # Calculate scaling factors
    scaling_factors = np.ceil(box_dimensions / unit_cell_lengths).astype(int)
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
right_supercell = create_rotated_supercell(unit_cell, right_box_dimensions, -right_rotation_angle) # negative means clockwise rotation

# Shift the right supercell to join with the left one
right_supercell.positions[:, 0] += left_box_dimensions[0]  # Shift in the x-direction

# Combine the two supercells
combined_supercell = left_supercell + right_supercell

# Trim any atoms that might be slightly outside the overall box due to rotations
final_positions = []
final_symbols = []

for atom in combined_supercell:
    pos = atom.position
    if (0 <= pos[0] < desired_box[0]) and (0 <= pos[1] < desired_box[1]) and (0 <= pos[2] < desired_box[2]):
        final_positions.append(pos)
        final_symbols.append(atom.symbol)

final_supercell = Atoms(positions=final_positions, symbols=final_symbols, cell=desired_box)
final_supercell.set_pbc(True)

# Visualize and save
view(final_supercell)
final_supercell.write('Silicon_grain_boundary_50x25x25.data', format='lammps-data')
final_supercell.write('Silicon_grain_boundary_50x25x25.xyz')