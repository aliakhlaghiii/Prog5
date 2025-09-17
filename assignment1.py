#!/usr/bin/env python3
"""Trapezoidal integration script.

This script approximates the integral of cos(x) over a given interval
using the trapezoidal rule, and compares the result with the exact
value of the integral.
"""

import math
import argparse


def trapezoid(func, lower, upper, steps=256):
    """Approximate the integral of a function using the trapezoidal rule.

    Args:
        func (Callable[[float], float]): Function to integrate.
        lower (float): Lower bound of integration.
        upper (float): Upper bound of integration.
        steps (int, optional): Number of subdivisions. Defaults to 256.

    Returns:
        float: Approximation of the integral.
    """
    step_size = (upper - lower) / steps
    integral = func(lower) + func(upper)

    for index in range(1, steps):
        integral += 2 * func(lower + index * step_size)

    return integral * (step_size / 2)


def main():
    """Parse command-line arguments and run trapezoidal integration."""
    parser = argparse.ArgumentParser(description="Trapezoidal integration demo")
    parser.add_argument("-a", "--lower", type=float, required=True,
                        help="Lower bound of integration")
    parser.add_argument("-b", "--upper", type=float, required=True,
                        help="Upper bound of integration")
    parser.add_argument("-n", "--steps", type=int, required=True,
                        help="Number of steps in the trapezoidal rule")
    args = parser.parse_args()

    func = math.cos
    exact = math.sin(args.upper) - math.sin(args.lower)
    approx = trapezoid(func, args.lower, args.upper, args.steps)
    error = abs(exact - approx)
    print(f"Steps: {args.steps}, Error: {error:.3f}")


if __name__ == "__main__":
    main()
