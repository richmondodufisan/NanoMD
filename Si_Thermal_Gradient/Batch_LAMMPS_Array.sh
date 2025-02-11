#!/bin/bash
#SBATCH --account=p32089  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=long  ### PARTITION (buyin, short, normal, etc)
#SBATCH --array=0-5         ## Adjust based on the number of LAMMPS scripts in the text file (0-indexed)
#SBATCH --nodes=1           ## Number of nodes per job
#SBATCH --ntasks-per-node=10 ## Number of CPUs per node
#SBATCH --time=168:00:00     ## Max runtime per job
#SBATCH --mem-per-cpu=1G    ## Memory per CPU
#SBATCH --job-name=Si_LAMMPS_Batch  ## Job name for identification
#SBATCH --output=log_%A_%a.lammps ## Output log file for each array task (%A is the job ID, %a is the array task ID)
#SBATCH --constraint="[quest8|quest9|quest10|quest11]"

# Load the necessary module
module purge
module load lammps/20200303-openmpi-4.0.5-intel-19.0.5.281

# Read the file containing LAMMPS script paths
IFS=$'\n' read -d '' -r -a scripts < SiliconSimulations.txt

# Get the current script to run based on the SLURM_ARRAY_TASK_ID
current_script=${scripts[$SLURM_ARRAY_TASK_ID]}

# Run the LAMMPS simulation
mpirun -np ${SLURM_NTASKS} lmp -in $current_script
