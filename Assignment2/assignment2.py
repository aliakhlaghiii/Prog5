#!/usr/bin/env python3
"""
Parallel numerical integration with the composite trapezoidal rule (mpi4py).

Usage (example):
  srun --ntasks=4 python3 assignment2.py -a 0 -b 1 -n 1000000
"""

import argparse as ap
import math
from mpi4py import MPI


# ---------- 1) The pure (serial) math core ----------
def f(x: float) -> float:
    """Function to integrate. You can hardcode it for the assignment."""
    return math.cos(x)


def trapezoid_serial(fn, a: float, b: float, n: int) -> float:
    """
    Composite trapezoidal rule on [a,b] with n subintervals.

    h = (b-a)/n
    Integral ≈ h * ( 0.5*f(a) + f(a+h) + f(a+2h) + ... + f(b-h) + 0.5*f(b) )
    """
    if n <= 0:
        raise ValueError("n must be positive")
    if a == b:
        return 0.0

    h = (b - a) / n
    s = 0.5 * fn(a) + 0.5 * fn(b)

    # add interior points: i = 1..n-1
    for i in range(1, n):
        x = a + i * h
        s += fn(x)

    return h * s


# ---------- 2) MPI wrapper around the math core ----------
def split_steps(n: int, size: int, rank: int):
    """
    Split n steps across 'size' ranks, handling remainder safely.

    Returns:
      start_step (inclusive), end_step (exclusive), local_n = end_step - start_step

    We split subinterval indices [0..n-1] into chunks.
    """
    base = n // size
    rem = n % size

    # ranks < rem get one extra step
    local_n = base + (1 if rank < rem else 0)

    # starting index depends on how many extra steps were given to earlier ranks
    start = rank * base + min(rank, rem)
    end = start + local_n
    return start, end, local_n


def trapezoid_parallel(fn, a: float, b: float, n: int, comm: MPI.Comm) -> float:
    """
    Parallel trapezoid: each rank integrates over its own sub-range.
    Then we sum partial integrals with Reduce.

    Key idea:
      Each rank gets a contiguous block of steps (subintervals).
      Convert that to a local [a_local, b_local] and local_n.
    """
    rank = comm.Get_rank()
    size = comm.Get_size()

    if n <= 0:
        raise ValueError("n must be positive")
    if a == b:
        return 0.0

    # global step size
    h = (b - a) / n

    # decide which subintervals this rank owns
    start_step, end_step, local_n = split_steps(n, size, rank)

    # map owned subinterval indices to a local [a_local, b_local]
    a_local = a + start_step * h
    b_local = a + end_step * h

    # compute local integral serially on local segment
    local_I = trapezoid_serial(fn, a_local, b_local, local_n)

    # sum all local integrals into rank 0
    total_I = comm.reduce(local_I, op=MPI.SUM, root=0)
    return total_I  # on non-root ranks, this is None


# ---------- 3) main: parse args, setup MPI, choose serial/parallel ----------
def main():
    parser = ap.ArgumentParser(description="MPI trapezoid integration")
    parser.add_argument("-a", type=float, required=True, help="Lower bound")
    parser.add_argument("-b", type=float, required=True, help="Upper bound")
    parser.add_argument("-n", type=int, required=True, help="Number of steps")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # (optional) ensure consistent timing start
    comm.Barrier()
    t0 = MPI.Wtime()

    if size == 1:
        # pure serial run
        result = trapezoid_serial(f, args.a, args.b, args.n)
    else:
        # parallel run
        result = trapezoid_parallel(f, args.a, args.b, args.n, comm)

    comm.Barrier()
    t1 = MPI.Wtime()

    # only rank 0 prints final result + timing
    if rank == 0:
        print(f"size={size} integral ≈ {result}")
        print(f"elapsed_seconds={t1 - t0}")


if __name__ == "__main__":
    main()
