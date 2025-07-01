#!/bin/bash

# This script shows how to use atomsk to build
# a supercell of carbon diamond

# First create the unit cell, then duplicate it,
# and output to XYZ and CFG for visualization
atomsk --create diamond 3.567 C -dup 4 4 1 C_diamond_supercell xsf cfg
