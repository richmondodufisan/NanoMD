#!/bin/bash
#SBATCH --account=p32089  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=short  ### PARTITION (buyin, short, normal, etc)
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=10 ## how many cpus or processors do you need on each computer
#SBATCH --time=01:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=2G ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=silicon_box  ## When you run squeue -u NETID this is how you can identify the job
#SBATCH --output=ase_log ## standard out and standard error goes to this file
#SBATCH --constraint="[quest8|quest9|quest10|quest11]"

mamba activate myenv

python3 silicon_box.py
