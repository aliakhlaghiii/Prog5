#!/usr/bin/python
import math
import sys
import argparse

def trapezoid(f, a, b, n=256):
    ''' function to implement integration based on lower bound (a),
    upper bound(b) and number of steps(n)'''
    
    h = (b - a) / n  
    I = (f(a) + f(b))  

    for index in range(1, n):
        I += 2 * f(a + index * h)  
    
    return I * (h / 2) 

def main():
    parser = argparse.ArgumentParser(description='Trapezoid model')
    parser.add_argument('-a', type=float, required=True, help='lower bound')
    parser.add_argument('-b', type=float, required=True, help='upper bound')
    parser.add_argument('-n', type=float, required=True, help='number of steps')
    args = parser.parse_args()
    
    f = math.cos
    exact_ans = math.sin(args.b) - math.sin(args.a)
    approx_ans = trapezoid(f, args.a, args.b, args.n)
    error = abs(exact_ans - approx_ans)
    print(f'{args.n}, {error:.3f}')
        