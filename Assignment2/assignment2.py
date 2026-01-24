#!/usr/bin/env python3
import argparse as ap
from mpi4py import MPI


# ---------------- 1) The pure (serial) core ----------------
def summation_serial(a: int, b: int) -> int:
    """Sum_{i=a}^{b-1} (i^2 + 1), b is exclusive."""
    s = 0
    for i in range(a, b):
        s += i * i + 1
    return s


# ---------------- 2) Split work across WORKERS ----------------
def split_steps(a: int, b: int, workers: int, worker_id: int):
    """
    Split [a,b) across 'workers' worker processes (NOT including rank0).
    worker_id is in [0 .. workers-1].
    Returns (start, end, local_n) with end exclusive.
    """
    n = b - a
    if n < 0:
        raise ValueError("b must be >= a")
    if workers <= 0:
        return a, a, 0  # no work

    base = n // workers
    rem = n % workers

    local_n = base + (1 if worker_id < rem else 0)
    start = a + worker_id * base + min(worker_id, rem)
    end = start + local_n
    return start, end, local_n


# ---------------- 3) MPI wrapper (parallel) ----------------
def summation_parallel(a: int, b: int, comm: MPI.Comm):
    """
    rank0 does no work, only reduces.
    workers (rank>=1) compute partial sums and Reduce(SUM) to rank0.
    """
    rank = comm.Get_rank()
    size = comm.Get_size()

    workers = size - 1

    if rank == 0:
        local_sum = 0
    else:
        worker_id = rank - 1  # map ranks 1..size-1 -> worker ids 0..workers-1
        start, end, _ = split_steps(a, b, workers, worker_id)
        local_sum = summation_serial(start, end)

    total_sum = comm.reduce(local_sum, op=MPI.SUM, root=0)
    return total_sum  # only rank0 gets the number; others get None


# ---------------- 4) main ----------------
def main():
    parser = ap.ArgumentParser(description="MPI summation of (i^2 + 1) over [a, b)")
    parser.add_argument("-a", type=int, required=True, help="start (inclusive)")
    parser.add_argument("-b", type=int, required=True, help="end (exclusive)")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    comm.Barrier()
    t0 = MPI.Wtime()

    if size == 1:
        result = summation_serial(args.a, args.b)
    else:
        result = summation_parallel(args.a, args.b, comm)

    comm.Barrier()
    t1 = MPI.Wtime()

    if rank == 0:
        print(f"size={size} sum={result}")
        print(f"elapsed_seconds={t1 - t0}")


if __name__ == "__main__":
    main()
