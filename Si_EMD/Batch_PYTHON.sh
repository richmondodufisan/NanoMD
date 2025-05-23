#!/bin/bash
#SBATCH --account=p32089  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=normal  ### PARTITION (buyin, short, normal, etc)
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=20 ## how many cpus or processors do you need on each computer
#SBATCH --time=24:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=100M ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=python  ## When you run squeue -u NETID this is how you can identify the job


module purge
module load mamba/24.3.0
module load git
eval "$(conda shell.bash hook)"

conda activate myenv

# python3 create_trajectory_video.py
python3 green_kubo_heatmap_HPC.py

# python3 green_kubo_calc.py
