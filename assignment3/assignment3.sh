#!/bin/bash
#SBATCH --job-name=assigment3
#SBATCH --partition=assemblix
#SBATCH --ntasks=33
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

STEPS=200000               
PY_SCRIPT=./assignment3.py
RESULTS=results.csv

if [[ ! -f "$RESULTS" ]]; then
  echo "workers,procs,steps,elapsed_s,user_s,sys_s,max_rss_kb,jobid,nodelist" > "$RESULTS"
fi

for workers in $(seq 1 32); do
  procs=$((workers + 1))
  /usr/bin/time -f "${workers},${procs},${STEPS},%e,%U,%S,%M,${SLURM_JOB_ID},${SLURM_JOB_NODELIST}" \
    -a -o "$RESULTS" \
    srun -n "$procs" /usr/bin/python3 "$PY_SCRIPT" -n "$STEPS"
done
