#!/usr/bin/python3
import math
import argparse
from mpi4py import MPI


def main():
    # mpi setup
    comm = MPI.COMM_WORLD   # gather all the ranks in a communicator 
    rank = comm.Get_rank()  # id of this rank
    size = comm.Get_size()  # total number of ranks (processes)

    # before broadcast, we have to define args on non-root ranks
    a = b = None
    n = None
    
    if rank == 0:
        # rank0 is in charge of control and broadcast 
        parser = argparse.ArgumentParser(description="Trapezoidal integration with MPI")
        parser.add_argument("-a", "--lower", type=float, required=True,
                            help="Lower bound of integration")
        parser.add_argument("-b", "--upper", type=float, required=True,
                            help="Upper bound of integration")
        parser.add_argument("-n", "--steps", type=int, required=True,
                            help="Number of steps in the trapezoidal rule")
        args = parser.parse_args()

        a, b, n = args.lower, args.upper, args.steps
        if n < 1:  # to provide a likely error for number of steps
            raise ValueError("steps (-n) must be >= 1")

    # broadcast arguments from rank0 to all ranks
    a, b, n = comm.bcast((a, b, n), root=0)

    # split the number of steps for each rank
    base = n // size  # minimum value of steps for each rank
    r = n % size   # residual -> given to the first r ranks
    local_steps = base + (1 if rank < r else 0)   # adding residuals to n first ranks
    i_start = rank * base + min(rank, r)   # index of start for rank
    i_end = i_start + local_steps   # index of end 

    # trapezoid rule
    # recall: trapezoid rule = (h/2)*f(a + 2*sigma f(a+ih) + f(b))
    h = (b - a) / n
    local_sum = 0.0   # initial each local rank

    for i in range(i_start + 1, i_end):   # without boundaries
        x = a + i * h   # f(x) = f(a + ih)
        local_sum += 2.0 * math.cos(x)
        
    # include left chunk-boundary (internal) with weight=2, only if we have work
    if local_steps > 0 and i_start > 0:
        local_sum += 2.0 * math.cos(a + i_start * h)

    # include global endpoints exactly once (only ranks touching them)
    if local_steps > 0 and i_start == 0:
        local_sum += math.cos(a)
    if local_steps > 0 and i_end == n:
        local_sum += math.cos(b)

    # call mpi.sum to sum(reduce) all the results from rank0
    total_sum = comm.reduce(local_sum, op=MPI.SUM, root=0)

    # rank0 has to report
    if rank == 0:
        approx = (h / 2.0) * total_sum
        exact = math.sin(b) - math.sin(a)
        error = abs(exact - approx)
        print("Controller (rank 0):")
        print(f"Ranks: {size}, Steps: {n}")
        print(f"Approx: {approx:.5f}, Exact: {exact:.5f}, Error: {error:.5e}")


if __name__ == "__main__":
    main()
