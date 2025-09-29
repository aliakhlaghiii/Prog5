# Assignment 3 — Profiling Python Code

## 1.1 Profiling results

**Q1. Can you see where the most calculation time is spent? What line number/function is it?**  
- Yes. The highest `tottime` is in the `sieve` function (line 2), which contains the nested loops that mark off composite numbers.  
- The second largest share is in the list comprehension `<listcomp>` (line 9), which scans the boolean array and collects the prime numbers.  
- `main` and `<module>` have very low `tottime` but high `cumtime`, because the entire program execution passes through them.

---

**Q2. What is the distribution of execution times in your file? Is it homogeneous (uniform) or is there a power distribution?**  
- The distribution is not uniform; it is right skewed.  
- Most of the execution time is concentrated in just two places; the computation inside `sieve` and the final list comprehension.  
- Input/output (`print`) and module imports take a bit time.
