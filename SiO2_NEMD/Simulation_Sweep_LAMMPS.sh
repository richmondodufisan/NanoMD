#!/bin/bash

# Base scaling factor settings
min_scale=1.0
max_scale=50.0
step=2.5

# File paths
base_python_script="silica_box.py"  # Replace with the Python script name
base_lammps_script="SiO2_Bulk_Hot_Cold.lammps"  # Replace with your LAMMPS script name
geometry_file_prefix="Silica_supercell"  # Prefix for generated geometry files
xyz_file_prefix="Silica_supercell"  # Prefix for XYZ files
output_dir="scaled_simulations"  # Directory to store results

# Create the output directory
mkdir -p $output_dir

# List to store the generated LAMMPS script paths
lammps_scripts=()

# Loop over scaling factors
for scale in $(seq $min_scale $step $max_scale); do
  # Scale dimensions
  scaled_x=$(echo "200 * $scale" | bc)
  scaled_y=25
  scaled_z=25
  
  # Generate scaled geometry
  scaled_python_script="${output_dir}/silicon_box_${scale}.py"
  scaled_geometry_file="${output_dir}/${geometry_file_prefix}_${scale}.data"
  scaled_xyz_file="${output_dir}/${xyz_file_prefix}_${scale}.xyz"
  
  # Copy and modify the Python script
  cp $base_python_script $scaled_python_script
  sed -i "s/desired_box = np.array(\[.*\])/desired_box = np.array(\\[$scaled_x, $scaled_y, $scaled_z\\])/" $scaled_python_script
  sed -i "s|trimmed_supercell\.write('.*\.data'|trimmed_supercell.write('${scaled_geometry_file}'|" $scaled_python_script
  sed -i "s|trimmed_supercell\.write('.*\.xyz'|trimmed_supercell.write('${scaled_xyz_file}'|" $scaled_python_script
  
  # Run the Python script to generate the geometry
  python3 $scaled_python_script
  
  # Adjust the LAMMPS input script
  scaled_lammps_script="${output_dir}/scaled_lammps_script_${scale}.lammps"
  cp $base_lammps_script $scaled_lammps_script

  sed -i "s/variable lx equal 200/variable lx equal $scaled_x/" $scaled_lammps_script
  sed -i "s/variable ly equal 25/variable ly equal $scaled_y/" $scaled_lammps_script
  sed -i "s/variable lz equal 25/variable lz equal $scaled_z/" $scaled_lammps_script

  # Update file paths, regions, and run time in the LAMMPS script
  sed -i "s|read_data .*|read_data ${scaled_geometry_file}|" $scaled_lammps_script
  
  # Adjust the run time based on the scale (longer for larger boxes)
  run_time=$(printf "%.0f" $(echo "1000000 * $scale" | bc))
  
  # Replace ONLY the LAST "run" command
  sed -i -E '$s/run [0-9]*/run '${run_time}'/' $scaled_lammps_script
  
  # Modify output file names to include scaling factor
  sed -i "s/temp_profile.txt/temp_profile_${scale}.txt/" $scaled_lammps_script
  
  sed -i "s/heat_flux_output_middle_left.txt/heat_flux_output_middle_left_${scale}.txt/" $scaled_lammps_script
  sed -i "s/heat_flux_output_middle_right.txt/heat_flux_output_middle_right_${scale}.txt/" $scaled_lammps_script

  sed -i "s/dump.lammpstrj/dump_${scale}.lammpstrj/" $scaled_lammps_script

  # Save the new LAMMPS script to the list
  lammps_scripts+=($scaled_lammps_script)
done

# Write the LAMMPS scripts to the text file, each on a new line
printf "%s\n" "${lammps_scripts[@]}" > SiliconSimulations.txt

# Update the SLURM job array range in the Batch_LAMMPS_Array.sh script
num_scripts=${#lammps_scripts[@]}
array_range="0-$(($num_scripts - 1))"

# Update the --array parameter in the Batch_LAMMPS_Array.sh file
sed -i "s/^#SBATCH --array=[^ ]*/#SBATCH --array=$array_range/" Batch_LAMMPS_Array.sh

echo "Simulations complete. Check the $output_dir directory for results."
echo "SLURM submission script Batch_LAMMPS_Array.sh has been updated with the correct job array range."