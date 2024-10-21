from ase.io import read, write

# Read CIF file
atoms = read("Si.cif")

# Convert to LAMMPS data file
write("Si.data", atoms, format="lammps-data")
