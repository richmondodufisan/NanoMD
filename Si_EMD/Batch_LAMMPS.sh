#!/bin/bash
#SBATCH --account=p32089  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=normal  ### PARTITION (buyin, short, normal, etc)
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=15 ## how many cpus or processors do you need on each computer
#SBATCH --time=48:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=200M ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=Si_EMD  ## When you run squeue -u NETID this is how you can identify the job
#SBATCH --exclude=qnode0565

module purge
module load lammps/20200303-openmpi-4.0.5-intel-19.0.5.281

mpirun -np ${SLURM_NTASKS} lmp -in Si_EMD.lammps

#srun --mpi=pmix_v3 lmp -in in.lj
