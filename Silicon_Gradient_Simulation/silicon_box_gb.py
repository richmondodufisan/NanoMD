import numpy as np
from ase.io import read
from ase.build import make_supercell
from ase import Atoms
from ase.visualize import view

## NEEDS WORK TO REFINE GB!!!!!!
## ROUGH CODE FOR NOW

# Load the unit cell from the CIF file
unit_cell = read('Si.cif')  # Replace with your CIF file path

# Desired box dimensions [x, y, z]
desired_box = np.array([50, 25, 25])

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

# Split the supercell into two parts for the grain boundary at x = 25
grain_boundary_position = 25.0

# Separate atoms into two regions based on the grain boundary
atoms_left = []
atoms_right = []
symbols_left = []
symbols_right = []

for atom in supercell:
    pos = atom.position
    if pos[0] < grain_boundary_position:
        atoms_left.append(pos)
        symbols_left.append(atom.symbol)
    else:
        atoms_right.append(pos)
        symbols_right.append(atom.symbol)

# Create Atoms object for the right half to apply rotation
atoms_right_obj = Atoms(positions=atoms_right, symbols=symbols_right, cell=desired_box)

# Rotate the right half (atoms_right_obj) by 60 degrees around the z-axis
rotation_angle = 60  # degrees
atoms_right_obj.rotate(rotation_angle, 'z', center=(grain_boundary_position, 0, 0))

# Combine atoms back together after rotation
all_positions = np.vstack([atoms_left, atoms_right_obj.get_positions()])
all_symbols = symbols_left + list(atoms_right_obj.get_chemical_symbols())

# Create a new Atoms object with the grain boundary structure
grain_boundary_supercell = Atoms(positions=all_positions, symbols=all_symbols, cell=desired_box)
grain_boundary_supercell.set_pbc(True)  # Enable periodic boundary conditions

# Visualize the structure to check the result
view(grain_boundary_supercell)

# Save the final structure to a LAMMPS data file
grain_boundary_supercell.write('Silicon_supercell_50x25x25_grain_boundary.data', format='lammps-data')

# Optionally, save as an XYZ file if needed
grain_boundary_supercell.write('Silicon_supercell_50x25x25_grain_boundary.xyz')
