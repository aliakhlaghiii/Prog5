# Assignment 3 — Profiling Python Code
## 1.1

**Profiling method:** I ran `cProfile` on `assignment3.py` and saved the output to `sieve_profiler_output.txt` using:
`python3 -m cProfile -s tottime assignment3.py > sieve_profiler_output.txt`

## 1.2 Profiling results

**Q1. Can you see where the most calculation time is spent? What line number/function is it?**  
- Yes. The highest `tottime` is in the `sieve` function (line 2), which contains the nested loops that mark off composite numbers.  
- The second largest share is in the list comprehension `<listcomp>` (line 9), which scans the boolean array and collects the prime numbers.  
- `main` and `<module>` have very low `tottime` but high `cumtime`, because the entire program execution passes through them.

---

**Q2. What is the distribution of execution times in your file? Is it homogeneous (uniform) or is there a power distribution?**  
- The distribution is not uniform; it is right skewed.  
- Most of the execution time is concentrated in just two places; the computation inside `sieve` and the final list comprehension.  
- Input/output (`print`) and module imports take a bit time.

1) At what point (cores) does it not matter anymore if you add any more (where does the plateau begin)?
Around 5 tasks. The best time is at 5 tasks; from 10–32 tasks the performance flattens or gets worse.
2) Can you explain this behavior?
As per-rank work decreases, MPI overheads (broadcast/gather/sync) and fixed serial parts dominate. The result is diminishing returns and a U-shaped runtime curve.
3) Does this behavior change when you run it on multiple hosts rather than just multiple cpu's on one host? Why or why not?
Yes, it gets worse. With --ntasks-per-node=1, adding nodes doesn’t reduce computation but adds inter-node latency/bandwidth costs, so network communication outweighs any benefit.