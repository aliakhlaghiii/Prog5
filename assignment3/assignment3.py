import math
def sieve(n=200000):
    A = [True] * (n+1) 
    A[0] = A[1] = False
    for i in range(2, math.isqrt(n)+1):
        if A[i]:
            for j in range(i*i, n+1, i):
                A[j] = False
    return [k for k in range(2, n+1) if A[k]]

       

def main():
    primes = sieve()
    print(primes)
    
if __name__ == '__main__':
    main()
    
    