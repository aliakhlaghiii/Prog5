#!/bin/bash
#SBATCH --job-name=assignment5
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=assignment5_%j.out
#SBATCH --error=assignment5_%j.err

python3 assignment5.py "$@"
