import matplotlib.pyplot as plt
from ase.io import read
from ase import Atoms
from ase.visualize.plot import plot_atoms

# Helper function to convert numeric types to element symbols
def load_xyz_with_types(filename):
    atoms = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        natoms = int(lines[0].strip())  # First line is the number of atoms
        positions = []
        types = []
        for line in lines[2:2 + natoms]:  # Skip first two lines
            atom_type, x, y, z = line.split()
            types.append(int(atom_type))  # Convert type to int
            positions.append([float(x), float(y), float(z)])
        
        # Map types to elements (e.g., 'H' for type 1, 'He' for type 2)
        type_to_symbol = {1: 'H', 2: 'He'}
        symbols = [type_to_symbol[t] for t in types]
        
        # Create ASE Atoms object
        atoms = Atoms(symbols=symbols, positions=positions)
    return atoms

# Load the modified XYZ data
atoms = load_xyz_with_types("initial_positions.xyz")

# Create a color map for atom types
atom_colors = {
    'H': "blue",  # Pseudo-element 'H' (Type 1)
    'He': "red",  # Pseudo-element 'He' (Type 2)
}

# Assign colors to atoms based on their types
colors = [atom_colors.get(atom.symbol, "gray") for atom in atoms]

# Visualize the atoms
fig, ax = plt.subplots()
plot_atoms(atoms, ax, radii=2.0, rotation=('90x,90y,0z'), colors=colors)

plt.axis("off")
plt.show()
