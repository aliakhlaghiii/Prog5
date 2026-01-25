#!/usr/bin/env python3
"""
Parallel numerical integration with the composite trapezoidal rule (mpi4py),
but with BROADCAST: rank 0 reads/owns the parameters and broadcasts them
to all ranks.

Usage:
  srun --ntasks=4 python3 assignment2.py -a 0 -b 1 -n 1000000
"""

import argparse as ap
import math
from mpi4py import MPI


# ---------- 1) The pure (serial) math core ----------
def f(x: float) -> float:
    return math.cos(x)


def trapezoid_serial(fn, a: float, b: float, n: int) -> float:
    if n <= 0:
        raise ValueError("n must be positive")
    if a == b:
        return 0.0

    h = (b - a) / n
    s = 0.5 * fn(a) + 0.5 * fn(b)
    for i in range(1, n):
        s += fn(a + i * h)
    return h * s


# ---------- 2) chunking ----------
def split_steps(n: int, size: int, rank: int):
    base = n // size
    rem = n % size
    local_n = base + (1 if rank < rem else 0)
    start = rank * base + min(rank, rem)
    end = start + local_n
    return start, end, local_n


def trapezoid_parallel(fn, a: float, b: float, n: int, comm: MPI.Comm):
    rank = comm.Get_rank()
    size = comm.Get_size()

    if n <= 0:
        raise ValueError("n must be positive")
    if a == b:
        return 0.0

    h = (b - a) / n

    start_step, end_step, local_n = split_steps(n, size, rank)

    a_local = a + start_step * h
    b_local = a + end_step * h

    local_I = trapezoid_serial(fn, a_local, b_local, local_n)

    total_I = comm.reduce(local_I, op=MPI.SUM, root=0)
    return total_I  # rank0 gets float; others get None


# ---------- 3) main with bcast ----------
def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Only rank0 parses args, then broadcasts (a,b,n) to everyone
    if rank == 0:
        parser = ap.ArgumentParser(description="MPI trapezoid integration (with bcast)")
        parser.add_argument("-a", type=float, required=True, help="Lower bound")
        parser.add_argument("-b", type=float, required=True, help="Upper bound")
        parser.add_argument("-n", type=int, required=True, help="Number of steps")
        args = parser.parse_args()
        params = (args.a, args.b, args.n)
    else:
        params = None

    a, b, n = comm.bcast(params, root=0)

    comm.Barrier()
    t0 = MPI.Wtime()

    if size == 1:
        result = trapezoid_serial(f, a, b, n)
    else:
        result = trapezoid_parallel(f, a, b, n, comm)

    comm.Barrier()
    t1 = MPI.Wtime()

    if rank == 0:
        print(f"size={size} integral ≈ {result}")
        print(f"elapsed_seconds={t1 - t0}")


if __name__ == "__main__":
    main()
