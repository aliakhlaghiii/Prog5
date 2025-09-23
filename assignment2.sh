#!/bin/bash
#SBATCH --job-name=assignment2
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --ntasks=33   # 32 ranks + 1 supervisor
#SBATCH --mem=1G

A=0
B=1.5707   # ~ pi/2
STEPS=1000

PY_SCRIPT=./assignment2.py

RESULTS=results.csv

if [[ ! -f "$RESULTS" ]]; then
  echo "workers,procs,steps,elapsed_s,user_s,sys_s,max_rss_kb,jobid,nodelist" > "$RESULTS"
fi
# %e: execution real time in sec
# %U: user CPU time
# %S: cpu time per each process
# %M: max RAM of launcher process
for workers in $(seq 1 32); do
  procs=$((workers + 1))  # +1 supervisor/controller

  /usr/bin/time -f "${workers},${procs},${STEPS},%e,%U,%S,%M,${SLURM_JOB_ID},${SLURM_JOB_NODELIST}" \
    -ao "$RESULTS" \
    srun --mpi=pmix -n "$procs" /usr/bin/python3 "$PY_SCRIPT" -a "$A" -b "$B" -n "$STEPS"
done

echo "Done. Results saved to ${RESULTS}"
