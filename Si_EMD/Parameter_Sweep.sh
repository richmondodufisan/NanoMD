#!/bin/bash

# Define parameter values
sample_vals=("20")
interval_vals=("50")

length_vals=("25" "50" "100" "200" "300")
temp_vals=("100" "300" "500" "700" "1000")

# Create the output directory if it doesn't exist
mkdir -p ./param_sweep

# Initialize the job list file
rm -f Si_EMD_List.txt
> Si_EMD_List.txt

# Counter for the array job index
job_index=0

# Loop over d first to generate the data files only once per d value
for d in "${length_vals[@]}"; do
    
    # Define data file names
    data_file="Silicon_supercell_${d}.data"
    xyz_file="Silicon_supercell_${d}.xyz"
    
    # Create a modified copy of silicon_box.py and run it to generate the required files
    cp silicon_box.py ./silicon_box_${d}.py
    sed -i "s/desired_box = np.array(\[109, 109, 109\])/desired_box = np.array([${d}, ${d}, ${d}])/" ./silicon_box_${d}.py
    sed -i "s|trimmed_supercell.write('Silicon_supercell.data'|trimmed_supercell.write('${data_file}'|" ./silicon_box_${d}.py
    sed -i "s|trimmed_supercell.write('Silicon_supercell.xyz'|trimmed_supercell.write('${xyz_file}'|" ./silicon_box_${d}.py
    python3 ./silicon_box_${d}.py
    
    # Loop over the remaining parameters
    for s in "${sample_vals[@]}"; do
        for p in "${interval_vals[@]}"; do
            for T in "${temp_vals[@]}"; do
                
                # Define output file names
                script_name="Si_EMD_sample_${s}_interval_${p}_length_${d}_temp_${T}.lammps"
                script_path="./param_sweep/$script_name"
                eq_log="equilibrium_sample_${s}_interval_${p}_length_${d}_temp_${T}.log"
                tc_log="thermal_conductivity_sample_${s}_interval_${p}_length_${d}_temp_${T}.log"
                heat_flux_file="heat_flux_sample_${s}_interval_${p}_length_${d}_temp_${T}.dat"
                j0jt_file="J0Jt_sample_${s}_interval_${p}_length_${d}_temp_${T}.dat"

                # Create a new LAMMPS script with modified parameters
                sed "s/^variable    T equal .*/variable    T equal ${T}/; \
                     s/^variable    p equal .*/variable    p equal ${p}/; \
                     s/^variable    s equal .*/variable    s equal ${s}/; \
                     s|log equilibrium.log|log ${eq_log}|; \
                     s|log thermal_conductivity.log|log ${tc_log}|; \
                     s|fix flux_output all print 1 \"\${step} \${Jx_raw} \${Jy_raw} \${Jz_raw}\" file heat_flux.dat screen no|fix flux_output all print 1 \"\${step} \${Jx_raw} \${Jy_raw} \${Jz_raw}\" file ${heat_flux_file} screen no|; \
                     s|fix JJ all ave/correlate \$s \$p \$d c_flux\[1\] c_flux\[2\] c_flux\[3\] type auto file J0Jt.dat ave running|fix JJ all ave/correlate \$s \$p \$d c_flux\[1\] c_flux\[2\] c_flux\[3\] type auto file ${j0jt_file} ave running|; \
                     s|read_data Silicon_supercell.data|read_data ${data_file}|" Si_EMD.lammps > "$script_path"
                
                # Append the script path to the job list file
                echo "$script_path" >> Si_EMD_List.txt
                
                # Increment job index
                ((job_index++))
            done
        done
    done
done

# Update the job array script with the correct array size
sed -i "s/^#SBATCH --array=0-[0-9]*/#SBATCH --array=0-$((job_index - 1))/" Batch_LAMMPS_Array.sh

echo "Generated $job_index LAMMPS input files and updated Batch_LAMMPS_Array.sh."
