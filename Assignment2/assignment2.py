#!/usr/bin/env python3
"""
Parallel numerical integration with the composite trapezoidal rule (mpi4py).
Rank 0 is a pure controller: it does NO integration work.

Usage:
  srun --ntasks=4 python3 assignment2.py -a 0 -b 1 -n 1000000
"""

import argparse as ap
import math
from mpi4py import MPI


# ---------- 1) The pure (serial) math core ----------
def f(x: float) -> float:
    """Hardcoded function to integrate."""
    return math.cos(x)


def trapezoid_serial(fn, a: float, b: float, n: int) -> float:
    """
    Composite trapezoidal rule on [a,b] with n subintervals.
    b is inclusive endpoint of the segment; n is the number of subintervals.
    """
    if n <= 0:
        return 0.0
    if a == b:
        return 0.0

    h = (b - a) / n
    s = 0.5 * fn(a) + 0.5 * fn(b)

    for i in range(1, n):
        x = a + i * h
        s += fn(x)

    return h * s


# ---------- 2) Split work across WORKERS (ranks 1..size-1) ----------
def split_steps_workers(n: int, workers: int, worker_id: int):
    """
    Split n steps across 'workers' worker processes (NOT including rank 0).
    worker_id is in [0 .. workers-1].

    Returns:
      start_step (inclusive), end_step (exclusive), local_n
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    if workers <= 0:
        return 0, 0, 0

    base = n // workers
    rem = n % workers

    local_n = base + (1 if worker_id < rem else 0)
    start = worker_id * base + min(worker_id, rem)
    end = start + local_n
    return start, end, local_n


def trapezoid_parallel_controller(fn, a: float, b: float, n: int, comm: MPI.Comm) -> float | None:
    """
    Rank 0: controller only (no computation)
    Ranks 1..size-1: compute local integrals
    Final result reduced to rank 0.
    """
    rank = comm.Get_rank()
    size = comm.Get_size()
    workers = size - 1

    if n <= 0 or a == b:
        # trivial integral
        local_I = 0.0
        total_I = comm.reduce(local_I, op=MPI.SUM, root=0)
        return total_I

    h = (b - a) / n

    if rank == 0:
        # controller does no work
        local_I = 0.0
    else:
        worker_id = rank - 1  # map ranks 1.. -> worker ids 0..
        start_step, end_step, local_n = split_steps_workers(n, workers, worker_id)

        # map step indices to the local [a_local, b_local]
        a_local = a + start_step * h
        b_local = a + end_step * h

        local_I = trapezoid_serial(fn, a_local, b_local, local_n)

    total_I = comm.reduce(local_I, op=MPI.SUM, root=0)
    return total_I  # only rank 0 gets the float; others get None


# ---------- 3) main ----------
def main():
    parser = ap.ArgumentParser(description="MPI trapezoid integration (rank0 controller only)")
    parser.add_argument("-a", type=float, required=True, help="Lower bound")
    parser.add_argument("-b", type=float, required=True, help="Upper bound")
    parser.add_argument("-n", type=int, required=True, help="Number of steps")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # timing (end-to-end)
    comm.Barrier()
    t0 = MPI.Wtime()

    if size == 1:
        # no workers exist; just run serial on the only process
        result = trapezoid_serial(f, args.a, args.b, args.n)
    else:
        result = trapezoid_parallel_controller(f, args.a, args.b, args.n, comm)

    comm.Barrier()
    t1 = MPI.Wtime()

    if rank == 0:
        print(f"size={size} integral ≈ {result}")
        print(f"elapsed_seconds={t1 - t0}")


if __name__ == "__main__":
    main()
