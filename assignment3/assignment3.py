#!/usr/bin/python3
"""Parallel segmented Sieve of Eratosthenes with mpi4py.

- Rank 0: parses args, computes base primes up to floor(sqrt(n)), broadcasts.
- Ranks 1..size-1: mark composites in their local [low..high] segment.
- Rank 0 gathers local prime lists and prints the count.
"""

import argparse
import math
from typing import List
from mpi4py import MPI

CONTROLLER_RANK = 0
DEBUG = False


def sieve(n: int) -> List[int]:
    """Return all primes in [2..n] using the classic sieve."""
    if n <= 1:
        raise ValueError('n must be > 1')
    mask = [True] * (n + 1)  # was A
    mask[0] = mask[1] = False
    for i in range(2, math.isqrt(n) + 1):
        if mask[i]:
            for j in range(i * i, n + 1, i):
                mask[j] = False
    return [k for k in range(2, n + 1) if mask[k]]


def main() -> None:
    """MPI orchestration: parse n, split work, gather results on rank 0."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == CONTROLLER_RANK:
        parser = argparse.ArgumentParser(
            description="Parallel sieve up to n"
        )
        parser.add_argument(
            "-n", "--limit", type=int, required=True,
            help="Upper bound n (find primes in [2..n])"
        )
        args = parser.parse_args()
        n_limit = args.limit  # was N
        if n_limit <= 1:
            raise ValueError("n must be > 1")
        if size < 2:
            raise RuntimeError(
                "Run with at least 2 ranks so rank 0 can orchestrate only."
            )
    else:
        n_limit = None  # type: ignore[assignment]

    # Broadcast N to everyone
    n_limit = comm.bcast(n_limit, root=CONTROLLER_RANK)

    if rank == CONTROLLER_RANK:
        base_primes = sieve(int(math.isqrt(n_limit)))
    else:
        base_primes = None
    base_primes = comm.bcast(base_primes, root=CONTROLLER_RANK)

    # use workers(rank0 is supervisor/controller)
    total = n_limit - 1                   # numbers from 2 to N
    workers = size - 1

    if rank == CONTROLLER_RANK or workers == 0:
        # rank 0 does no local compute
        local_primes: List[int] = []
        # (optional) debug
        # print(f"Rank 0: orchestrator only")
    else:
        w = rank - 1    # worker index in [0..workers-1]
        base, r = divmod(total, workers)    # base and residuals
        local_n = base + (1 if w < r else 0)
        offset = w * base + min(w, r)   # index of start in each rank

        low = 2 + offset   # start number of interval in each rank
        high = low + local_n - 1 if local_n > 0 else None   # end of interval

        if DEBUG:
            print(f"Rank {rank}: low={low}, high={high}, count={local_n}")

        if local_n > 0 and high is not None:
            seg = [True] * local_n   # making an array of Trues for each rank
            for p in base_primes:  # type: ignore[union-attr]
                # first multiple of p within [low..high], but not before p*p
                start = max(p * p, math.ceil(low / p) * p)
                if start > high:
                    continue
                for m in range(start, high + 1, p):
                    seg[m - low] = False   # mark multiples as composite
            local_primes = [i for i, ok in enumerate(seg, start=low) if ok]
        else:
            local_primes = []

    # Gather results to rank 0
    all_primes = comm.gather(local_primes, root=CONTROLLER_RANK)

    if rank == CONTROLLER_RANK:
        primes: List[int] = []
        for chunk in all_primes:
            primes.extend(chunk)
        print(f"number of primes = {len(primes)}")
        # print(primes)   # optional


if __name__ == '__main__':
    main()
