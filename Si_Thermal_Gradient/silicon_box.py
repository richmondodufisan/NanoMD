import numpy as np
from ase.io import read
from ase.build import make_supercell
from ase import Atoms
from ase.visualize import view

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

# Trim atoms that fall outside the desired box dimensions
trimmed_positions = []
trimmed_symbols = []

for atom in supercell:
    pos = atom.position
    if (0 < pos[0] < desired_box[0]) and (0 < pos[1] < desired_box[1]) and (0 < pos[2] < desired_box[2]):
        trimmed_positions.append(pos)
        trimmed_symbols.append(atom.symbol)

# Create a new Atoms object with only the atoms inside the box
trimmed_supercell = Atoms(positions=trimmed_positions, symbols=trimmed_symbols, cell=desired_box)
trimmed_supercell.set_pbc(True)  # Enable periodic boundary conditions

# Visualize the structure to check the result
# view(trimmed_supercell)

# Save the final trimmed structure to a LAMMPS data file
trimmed_supercell.write('Silicon_supercell_50x25x25.data', format='lammps-data')

# Optionally, you can also save as an XYZ file if needed
trimmed_supercell.write('Silicon_supercell_50x25x25.xyz')
