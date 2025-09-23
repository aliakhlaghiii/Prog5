#!/bin/bash
#SBATCH --job-name=assignment2
#SBATCH --partition=workstations
#SBATCH --ntasks=33
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

A=0
B=1.5707
STEPS=1000
PY_SCRIPT=./assignment2.py
RESULTS=results.csv

# if results.csv does not exist, create it and write the header line
if [[ ! -f "$RESULTS" ]]; then
  echo "workers,procs,steps,elapsed_s,user_s,sys_s,max_rss_kb,jobid,nodelist" > "$RESULTS"
fi

"""
use time to show the performance:
parameters:
    %e: execution real time(sec)
    %U: total CPU time spent in user mode (seconds)
    %S: total CPU time spent in kernel mode (seconds)
    %M: maximum resident set size (KB)
"""
for workers in $(seq 1 32); do
  procs=$((workers + 1))  
  /usr/bin/time -f "${workers},${procs},${STEPS},%e,%U,%S,%M,${SLURM_JOB_ID},${SLURM_JOB_NODELIST}" \
    -ao "$RESULTS" \   # -ao to append outputs one after another
    srun -n "$procs" /usr/bin/python3 "$PY_SCRIPT" -a "$A" -b "$B" -n "$STEPS"
done
