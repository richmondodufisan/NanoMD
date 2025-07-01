#!/bin/bash

# ======================
# CONFIGURATION SECTION
# ======================

CIF_FILE="Si_cubic.cif"         # Input CIF file (make sure it's in the same folder)
MISO_ANGLE="$1"           # Misorientation angle in degrees (e.g., 15)
DUP_X=20                  # ~5 Å lattice → 20 x ≈ 100 Å per grain
DUP_Y=5                   # ~5 x 5 = 25 Å → total ~50 Å
DUP_Z=5                   # same for Z

# ======================
# SCRIPT LOGIC SECTION
# ======================

if [[ ! -f "$CIF_FILE" ]]; then
    echo "Error: $CIF_FILE not found."
    exit 1
fi

if [[ -z "$MISO_ANGLE" ]]; then
    echo "Usage: bash make_silicon_gb.sh <misorientation_angle_in_degrees>"
    exit 1
fi

echo "➡ Creating left grain (0° rotation)..."
atomsk "$CIF_FILE" -rotate z 0 -duplicate $DUP_X $DUP_Y $DUP_Z left.xsf

echo "➡ Creating right grain (+$MISO_ANGLE° rotation)..."
atomsk "$CIF_FILE" -rotate z $MISO_ANGLE -duplicate $DUP_X $DUP_Y $DUP_Z right.xsf

echo "➡ Merging left and right grains..."
atomsk --merge X 2 left.xsf right.xsf Si_GB_200x50x50.lmp

echo "➡ Converting to XYZ for visualization..."
atomsk Si_GB_200x50x50.lmp Si_GB_200x50x50.xyz

echo "✅ Done! Output files:"
echo "   • Si_GB_200x50x50.lmp (LAMMPS)"
echo "   • Si_GB_200x50x50.xyz (XYZ for Ovito)"

