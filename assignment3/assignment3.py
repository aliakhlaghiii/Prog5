import math
def main():
    n = int(input('Enter a number:'))
    if n <= 1:
        raise ValueError('n must be > 1')
    A = [True] * (n+1) 
    A[0] = A[1] = False
    for i in range(2, math.isqrt(n)+1):
        if A[i]:
            for j in range(i*i, n+1, i):
                A[j] = False
    primes = [k for k in range(2, n+1) if A[k]]

    return primes        

if __name__ == '__main__':
    print(main())
    
    