"""Checks for 'The sharp small-noise constant for weighted majority'.

Run with Python 3 (standard library only). Exact rational arithmetic is used
where the paper states an identity or an exact value.
"""
from fractions import Fraction as Q
from itertools import product
from math import acos, comb, pi, sqrt
import random


def sgnp(u):
    return 1 if u >= 0 else -1


def sgn0(u):
    return (u > 0) - (u < 0)


def disagreement(f, n, eps):
    """Exact P(f(X) != f(N_eps X)) for f on {-1,1}^n (any number of values)."""
    pts = list(product((-1, 1), repeat=n))
    total = Q(0)
    for x in pts:
        fx = f(x)
        for flips in product((0, 1), repeat=n):
            y = tuple(-xi if fl else xi for xi, fl in zip(x, flips))
            if f(y) != fx:
                m = sum(flips)
                total += eps ** m * (1 - eps) ** (n - m)
    return total / len(pts)


def block_side(w, b, eps):
    """k^{-1} sum_j E[f(S) sgn(Z_j)] computed over all partitions (Lemma 2.1)."""
    n = len(w)
    k = int(1 / eps)
    probs = [1 - k * eps] + [eps] * k
    total = Q(0)
    for blocks in product(range(k + 1), repeat=n):
        pb = Q(1)
        for j in blocks:
            pb *= probs[j]
        if pb == 0:
            continue
        acc = Q(0)
        for x in product((-1, 1), repeat=n):
            s = sum(wi * xi for wi, xi in zip(w, x))
            fs = sgnp(s - b)
            for j in range(1, k + 1):
                z = sum(wi * xi for wi, xi, bj in zip(w, x, blocks) if bj == j)
                acc += fs * sgn0(z)
        total += pb * acc / 2 ** n
    return total / k


print("Lemma 2.1: N_eps(f) = k^{-1} sum_j E[f(S) sgn Z_j], exact")
cases = [((1, 1), 0, Q(1, 3)), ((1, 2, 3), 1, Q(2, 7)), ((1, -1, 3, 2), 1, Q(2, 5)),
         ((10, 1, 1), -10, Q(1, 4)), ((4, 0, -2), 2, Q(3, 8)), ((2, 2, 1), 0, Q(1, 3))]
for w, b, eps in cases:
    f = lambda x, w=w, b=b: sgnp(sum(wi * xi for wi, xi in zip(w, x)) - b)
    lhs = disagreement(f, len(w), eps)
    rhs = block_side(w, b, eps)
    assert lhs == rhs, (w, b, eps, lhs, rhs)
    print(f"  w={w} b={b} eps={eps}: {lhs} (agree)")


def abs_mean_thinned(w, eps):
    """E|V| for V = sum w_i X_i I_i, exact enumeration (floats)."""
    n = len(w)
    tot = 0.0
    for sel in product((0, 1), repeat=n):
        m = sum(sel)
        p = eps ** m * (1 - eps) ** (n - m)
        idx = [i for i in range(n) if sel[i]]
        acc = 0.0
        for x in product((-1, 1), repeat=len(idx)):
            acc += abs(sum(w[i] * xi for i, xi in zip(idx, x)))
        tot += p * acc / 2 ** len(idx)
    return tot


print("Lemma 3.3: E|W| <= sqrt(2/pi) (1-2 eps)^(-1/2), random weights")
rng = random.Random(20261119)
worst = 0.0
for trial in range(60):
    n = rng.randint(1, 10)
    w = [rng.gauss(0, 1) if rng.random() < 0.7 else rng.uniform(-3, 3) for _ in range(n)]
    nrm = sqrt(sum(v * v for v in w))
    w = [v / nrm for v in w]
    eps = rng.choice([0.02, 0.05, 0.1, 0.2, 0.3, 0.45])
    ew = abs_mean_thinned(w, eps) / sqrt(eps)
    bound = sqrt(2 / pi) / sqrt(1 - 2 * eps)
    assert ew <= bound + 1e-12, (w, eps, ew, bound)
    worst = max(worst, ew / bound)
print(f"  60 cases, largest ratio E|W| / bound = {worst:.6f}")


def maj_disagreement(n, eps):
    """N_eps(Maj_n) for odd n, by the distribution of the number of plus signs."""
    half = 2.0 ** (-n)
    tot = 0.0
    for a in range(n + 1):
        pa = comb(n, a) * half
        fx = 1 if 2 * a > n else -1
        pu = [comb(a, i) * eps ** i * (1 - eps) ** (a - i) for i in range(a + 1)]
        pv = [comb(n - a, j) * eps ** j * (1 - eps) ** (n - a - j) for j in range(n - a + 1)]
        ch = 0.0
        for i, p1 in enumerate(pu):
            for j, p2 in enumerate(pv):
                plus = a - i + j
                if (1 if 2 * plus > n else -1) != fx:
                    ch += p1 * p2
        tot += pa * ch
    return tot


print("Simple majority: N_eps / sqrt(eps) against 2/pi = %.9f" % (2 / pi))
print("     n      eps    N/sqrt(eps)   arccos(1-2eps)/(pi sqrt(eps))")
for n, eps in [(41, 1 / 41), (151, 1 / 151), (401, 1 / 401), (151, 1 / 25), (401, 1 / 25),
               (401, 1 / 100), (801, 1 / 100)]:
    r = maj_disagreement(n, eps) / sqrt(eps)
    lim = acos(1 - 2 * eps) / (pi * sqrt(eps))
    print(f"  {n:5d}  {eps:7.5f}   {r:.9f}   {lim:.9f}")

def maj_abs_v(n, eps):
    """E|V| for V = n^{-1/2} sum X_i I_i: the thinned count m ~ Bin(n, eps), then |2B_m - m|."""
    tot = 0.0
    for m in range(n + 1):
        pm = comb(n, m) * eps ** m * (1 - eps) ** (n - m)
        if pm < 1e-300:
            continue
        e = sum(comb(m, j) * abs(2 * j - m) for j in range(m + 1)) / 2.0 ** m
        tot += pm * e
    return tot / sqrt(n)


print("Proposition 3.2 for simple majority (b = 0): N/sqrt(eps) against 2 phi(0) E|V|/sqrt(eps)")
print("     n      eps    N/sqrt(eps)   main term    difference")
phi0 = 1 / sqrt(2 * pi)
prev = None
for n, eps in [(101, 1 / 25), (401, 1 / 25), (801, 1 / 25), (401, 1 / 100), (801, 1 / 100)]:
    r = maj_disagreement(n, eps) / sqrt(eps)
    mt = 2 * phi0 * maj_abs_v(n, eps) / sqrt(eps)
    d = abs(r - mt)
    print(f"  {n:5d}  {eps:7.5f}   {r:.6f}     {mt:.6f}    {d:.6f}")
    assert d < 0.1, (n, eps, r, mt)
    assert mt <= sqrt(2 / pi) * 2 * phi0 / sqrt(1 - 2 * eps) + 1e-12

print("Proposition 5.1 at n = 256, eps = 1/16")
n = 256
lb = (1 - Q(1, n * n)) ** n * (1 - (1 - Q(1, 16)) ** n)
assert lb > Q(99, 100) and lb > Q(996, 1000)
print(f"  lower bound = {float(lb):.9f} > 0.996; 2 sqrt(eps) = 1/2")
print("ALL CHECKS PASS")
