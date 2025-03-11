#!/bin/bash
#SBATCH --account=p32089  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=normal  ### PARTITION (buyin, short, normal, etc)
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=20 ## how many cpus or processors do you need on each computer
#SBATCH --time=12:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=1G ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=Si_thermal_grad  ## When you run squeue -u NETID this is how you can identify the job
#SBATCH --constraint="[quest8|quest9|quest10|quest11]"
#SBATCH --exclude=qnode0565,qnode0626,qnode0637,qnode0019,qnode0115,qnode1201

module purge
module load lammps/20200303-openmpi-4.0.5-intel-19.0.5.281

mpirun -np ${SLURM_NTASKS} lmp -in Si_EMD.lammps

#srun --mpi=pmix_v3 lmp -in in.lj
