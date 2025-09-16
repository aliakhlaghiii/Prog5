#!/bin/bash
#SBATCH --job-name=assignment1
#SBATCH --comment='implement trapezoid integration'

#SBATCH --account=Ali
#SBATCH --partition=workstations

#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

#SBATCH --mail-user=a.akhlaghi@st.hanze.nl
#SBATCH --mail-type=ALL

#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

A=0
B=1.5707   # ~ pi/2(radian)

# Create/reset results file with header
echo "n,error" > results.csv

# List of n values to test
Ns=(500 750 1000)

for n in "${Ns[@]}"; 
do
  # Run .py script on the allocated node
  line=$(srun /usr/bin/python3 ./assignment1.py -a "$A" -b "$B" -n "$n")
  echo "$line" >> results.csv
done

echo "Done. Results saved to results.csv)"
