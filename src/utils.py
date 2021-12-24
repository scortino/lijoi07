from functools import lru_cache
from scipy.special import comb, factorial

# TODO: these are trivial implementations straight from the paper, try to vectorize or speed up using numba

# rising factorial
def rising_factorial(a, n):
    p = a
    for j in range(1, n):
        p *= (a + j)
    return p

# Central generalized factorial coefficient
def generalized_factorial(n, k, sigma):
    if (n == k == 0):
        return 1
    s = 0
    for j in range(k + 1):
        s += (-1)**j * comb(k, j) * rising_factorial(-j * sigma, n)
    return s / factorial(k)

# Non-central generalized factorial coefficient
def noncentral_generalized_factorial(n, k, sigma, gamma):
    s = 0
    for j in range(k, n + 1):
        s += comb(n, j) * generalized_factorial(j, k, sigma) * rising_factorial(-gamma, n - j)
    return s

# Unsigned stirling number of the first kind
# https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
# n>=0, k>=0
@lru_cache
def stirling(n, k):
    if n == k == 0:
        return 1
    elif (n == 0) or (k == 0):
        return 0
    return (n - 1) * stirling(n - 1, k) + stirling(n - 1, k - 1)
