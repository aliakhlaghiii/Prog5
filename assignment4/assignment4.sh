#!/bin/bash
#SBATCH --job-name=assignment4
#SBATCH --partition=workstations
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

/usr/bin/time -p srun -u python3 assignment4.py
