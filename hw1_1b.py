from collections import defaultdict
from fractions import Fraction

def propensities(x1, x2, x3):
    # a1 = k1 * C(x1,2) * x2, k1=1
    a1 = Fraction(1, 2) * x1 * (x1 - 1) * x2 if (x1 >= 2 and x2 >= 1) else Fraction(0, 1)

    # a2 = k2 * x1 * C(x3,2), k2=2 => x1*x3*(x3-1)
    a2 = Fraction(x1 * x3 * (x3 - 1), 1) if (x1 >= 1 and x3 >= 2) else Fraction(0, 1)

    # a3 = k3 * x2 * x3, k3=3
    a3 = Fraction(3 * x2 * x3, 1) if (x2 >= 1 and x3 >= 1) else Fraction(0, 1)

    return a1, a2, a3

def step_distribution(dist):
    """
    dist: dict mapping (x1,x2,x3) -> probability (Fraction)
    returns next-step distribution
    """
    nd = defaultdict(Fraction)
    for (x1, x2, x3), p in dist.items():
        a1, a2, a3 = propensities(x1, x2, x3)
        a0 = a1 + a2 + a3

        if a0 == 0:
            nd[(x1, x2, x3)] += p  # stuck state
            continue

        if a1 > 0:
            nd[(x1 - 2, x2 - 1, x3 + 4)] += p * a1 / a0
        if a2 > 0:
            nd[(x1 - 1, x2 + 3, x3 - 2)] += p * a2 / a0
        if a3 > 0:
            nd[(x1 + 2, x2 - 1, x3 - 1)] += p * a3 / a0

    return nd

def mean_and_variance(dist, idx):
    """
    idx=0 for X1, 1 for X2, 2 for X3
    """
    EX  = sum(Fraction(state[idx], 1) * p for state, p in dist.items())
    EX2 = sum(Fraction(state[idx]**2, 1) * p for state, p in dist.items())
    Var = EX2 - EX * EX
    return EX, Var

def main():
    # initial distribution: state [9,8,7] with prob 1
    dist = {(9, 8, 7): Fraction(1, 1)}

    # propagate 7 steps
    for _ in range(7):
        dist = step_distribution(dist)

    m1, v1 = mean_and_variance(dist, 0)
    m2, v2 = mean_and_variance(dist, 1)
    m3, v3 = mean_and_variance(dist, 2)

    # print as decimals
    print(f"E[X1] = {float(m1):.6f}, Var(X1) = {float(v1):.6f}")
    print(f"E[X2] = {float(m2):.6f}, Var(X2) = {float(v2):.6f}")
    print(f"E[X3] = {float(m3):.6f}, Var(X3) = {float(v3):.6f}")

if __name__ == "__main__":
    main()