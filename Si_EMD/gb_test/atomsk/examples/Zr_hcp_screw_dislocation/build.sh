#!/bin/bash

rm -f Zr_screw.*

# Write elastic constants of hcp Zr into text file
echo "# Elastic constants of hcp zirconium" > Zr_Cij.txt
echo "# E. S. Fisher and C. J. Renken, Phys. Rev. 135, A482 (1964)." >> Zr_Cij.txt
echo "elastic Voigt" >> Zr_Cij.txt
echo "155.4 155.4 172.5  #C11 C22 C33" >> Zr_Cij.txt
echo " 64.6  64.6  68.2  #C23 C31 C12" >> Zr_Cij.txt
echo " 36.3  36.3  44.1  #C44 C55 C66" >> Zr_Cij.txt

# Create oriented unit cell of hcp Zr
# Read elastic constants from file with option "-properties"
# Duplicate it to form a 40x40x1 supercell
# Introduce a 1/3[2-1-10] screw dislocation in the center of the box
# Select atoms out of central cylinder of radius 60 A, and fix them
# Output final result to Zr_scrw.cfg for visualization (Atomeye or OVITO),
# and Zr_screw.lmp for simulation with LAMMPS
atomsk --create hcp 3.232 5.180896 Zr orient [0001] [0-110] [2-1-10] \
       -properties Zr_Cij.txt \
       -duplicate 40 40 1 \
       -dislocation 0.501*box 0.501*box screw Z Y 3.232  \
       -select out cylinder Z 0.5*box 0.5*box 60 -fix all \
       Zr_screw.cfg lmp -vvv

