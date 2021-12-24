from scipy.special import comb, factorial

# TODO: this are trivial implementations straight from the paper, try to vectorize or speed up using numba

# Ascending factorial
def ascending_factorial(a, n):
    p = a
    for j in range(1, n):
        p *= (a + j)
    return p

# Central generalized factorial coefficient
def generalized_factorial(n, k, sigma):
    if (n == 0) and (k == 0):
        return 1.0
    s = 0
    for j in range(k + 1):
        s += (-1)**j * comb(k, j) * ascending_factorial(-j * sigma, n)
    return s / factorial(k)

# Non-central generalized factorial coefficient
def noncentral_generalized_factorial(n, k, sigma, gamma):
    s = 0
    for j in range(k, n + 1):
        s += comb(n, j) * generalized_factorial(j, k, sigma) * ascending_factorial(-gamma, n - j)
    return s