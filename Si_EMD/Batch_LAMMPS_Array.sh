#!/bin/bash
#SBATCH --account=p32089  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=long  ### PARTITION (buyin, short, normal, etc)
#SBATCH --array=0-15
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15 ## how many cpus or processors do you need on each computer
#SBATCH --time=168:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=200M ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=emd_array  ## When you run squeue -u NETID this is how you can identify the job
#SBATCH --constraint="[quest8|quest9|quest10|quest11]"


module purge
module load lammps/20200303-openmpi-4.0.5-intel-19.0.5.281

IFS=$'\n' read -d '' -r -a lines < Si_EMD_List.txt

mpirun -np ${SLURM_NTASKS} lmp -in ${lines[$SLURM_ARRAY_TASK_ID]}

