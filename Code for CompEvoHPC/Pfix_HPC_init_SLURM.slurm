#!/bin/bash

#SBATCH --job-name=Pfix_Estimation
#SBATCH --ntasks=20
#SBATCH --nodes=1           
#SBATCH --mem-per-cpu=6gb                   
#SBATCH --time=20:00:00   
#SBATCH --partition=windfall
#SBATCH --account=windfall  

module load python/3.8
python3 HPC_Pfix_Simulation.py
