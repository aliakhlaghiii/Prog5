#!/bin/bash
#SBATCH --job-name=assignment6
#SBATCH --partition=assemblix
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

srun -u python3 assignment6.py