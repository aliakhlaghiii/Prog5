#!/usr/bin/python3
import math
import argparse
from mpi4py import MPI


def main():
    # --- MPI setup ---
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # --- rank0: parse args ---
    if rank == 0:
        parser = argparse.ArgumentParser(description="Parallel prime counting with MPI")
        parser.add_argument("-n", "--limit", type=int, required=True,
                            help="Upper bound n (find primes in [2..n])")
        args = parser.parse_args()
        limit = args.limit
        if limit < 2:
            raise ValueError("limit (-n) must be >= 2")

        # build base primes up to floor(sqrt(limit)) inline (classic sieve)
        lim = math.isqrt(limit)
        baseA = [True] * (lim + 1)
        if lim >= 0: baseA[0] = False
        if lim >= 1: baseA[1] = False
        for p in range(2, math.isqrt(lim) + 1):
            if baseA[p]:
                start = p * p
                step = p
                baseA[start:lim + 1:step] = [False] * (((lim - start) // step) + 1)
        base_primes = [i for i in range(2, lim + 1) if baseA[i]]
    else:
        limit = None
        base_primes = None

    # --- broadcast shared params ---
    limit = comm.bcast(limit if rank == 0 else None, root=0)
    base_primes = comm.bcast(base_primes if rank == 0 else None, root=0)

    # --- split work among ranks (rank0 = controller) ---
    workers = size - 1
    if workers <= 0:
        if rank == 0:
            print("need at least 2 ranks")
        return

    lower = 2
    N = limit - lower + 1  # total numbers in [2..limit]

    if rank == 0:
        local_list = None
        local_count = 0  # controller doesn't compute
    else:
        w = rank - 1                       # worker index in [0..workers-1]
        q, r = divmod(N, workers)          # golden formula
        len_w  = q + (1 if w < r else 0)
        start_w = lower + w * q + min(w, r)
        end_w   = start_w + len_w - 1      # inclusive

        if len_w <= 0:
            local_list = []
            local_count = 0
        else:
            # --- segmented sieve on [start_w .. end_w] ---
            low  = max(start_w, 2)
            high = end_w
            seg_len = high - low + 1
            mark = [True] * seg_len
            for p in base_primes:
                p2 = p * p
                # first multiple of p inside [low..high], never before p^2
                start = max(p2, ((low + p - 1) // p) * p)
                if start > high:
                    continue
                for j in range(start, high + 1, p):
                    mark[j - low] = False
            local_list = [low + i for i, ok in enumerate(mark) if ok]
            local_count = len(local_list)

    # --- reduce only: total count (no gather) ---
    total_count = comm.reduce(local_count, op=MPI.SUM, root=0)

    # ======= Added: show first/last K primes without heavy gather =======
    K = 10  # می‌تونی 20 هم بذاری
    local_head = [] if rank == 0 else local_list[:K]
    local_tail = [] if rank == 0 else local_list[-K:]

    heads = comm.gather(local_head, root=0)
    tails = comm.gather(local_tail, root=0)

    # --- rank0: report ---
    if rank == 0:
        # build first K by scanning heads in rank order
        first_k = []
        for rnk in range(1, size):
            for x in heads[rnk]:
                if len(first_k) < K:
                    first_k.append(x)
                else:
                    break
            if len(first_k) >= K:
                break

        # build last K by concatenating all tails then taking the last K
        all_tails = []
        for rnk in range(1, size):
            all_tails.extend(tails[rnk])
        last_k = all_tails[-K:] if len(all_tails) >= K else all_tails

        print(f"Total primes up to {limit}: {total_count} (via MPI.REDUCE)")
        print(f"First {min(K, total_count)} primes: {first_k}")
        print(f"Last  {min(K, total_count)} primes:  {last_k}")


if __name__ == "__main__":
    main()
