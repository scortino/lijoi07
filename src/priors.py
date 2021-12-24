import numpy as np
from scipy.special import comb
from decimal import Decimal
import sys

from .utils import rising_factorial, generalized_factorial, noncentral_generalized_factorial, stirling

# TODO: these are trivial implementations straight from the paper, try to vectorize or speed up using numba

class GibbsType:
    def __init__(self, sigma, V):
        self.sigma = sigma
        self.V = V
        self.C = lambda n, k: generalized_factorial(n, k, sigma)
        self.nC = lambda n, k, gamma: noncentral_generalized_factorial(n, k, sigma, gamma)

    def prior_prob_K(self, n, k):
        return self.V(n, k) * self.C(n, k) / (self.sigma**k)

    def posterior_prob_K(self, m, k, n, j):
        return self.V(n + m, j + k) * self.nC(m, k, -n + j * self.sigma) / (self.V(n, j) * (self.sigma**k))

    def Dhat_m(self, m, n, j):
        s = 0
        for k in range(m + 1):
            s += self.V(n + m + 1, j + k + 1) * self.nC(m, k, -n + j * self.sigma) / (self.sigma**k)
        return s / self.V(n, j)

class DirichletProcess(GibbsType):
    def __init__(self, theta):
        self.theta = theta
        super().__init__(0, lambda n, k: (theta**k) / rising_factorial(theta, n))

    def prior_prob_K(self, n, k):
        return self.theta**k * stirling(n, k) / rising_factorial(self.theta, n)

    # Notice j unused
    def posterior_prob_K(self, m, k, n):
        s = 0
        for l in range(k, m + 1):
            s += comb(m, l) * stirling(l, k) * rising_factorial(n, m - l)
        return s * self.theta**k * rising_factorial(self.theta, n) / rising_factorial(self.theta, n + m)
    
    # Notice j unused
    def Dhat_m(self, m, n):
        s1 = 0
        for k in range(m + 1):
            s2 = 0 
            for l in range(k, m + 1):
                s2 += comb(m, l) * stirling(l, k) * rising_factorial(n, m - l)
            s1 += self.theta**k * s2
        return s1 * self.theta / rising_factorial(self.theta + n, m + 1)

class PoissonDirichletProcess(GibbsType):
    def __init__(self, sigma, theta):
        self.theta = theta
        # TODO: find V expression for PD process
        super().__init__(sigma, None)
    
    def prior_prob_K(self, n, k):
        p = 1
        for i in range(1, k):
            p *= (self.theta + i * self.sigma)
            # TODO Expression is correct, but there is a problem of floating point precision when sigma is small and k is large
        return p * self.C(n, k) / (self.sigma)**k  / rising_factorial(self.theta + 1, n - 1)

    def posterior_prob_K(self, m, k, n, j):
        p = 1
        for i in range(j, j + k):
            p *= self.theta + i * self.sigma * self.nC(m, k, -n + j * self.sigma)
        return p * rising_factorial(self.theta + 1, n - 1) / ((self.theta)**k * rising_factorial(self.theta + 1, n + m))

    def Dhat_m(self, m, k, n, j):
        s = 0
        for k in range(m + 1):
            p = 1
            for i in range(j, j + k + 1):
                p *= self.theta + i * self.sigma * self.nC(m, k, -n + j * self.sigma)
            s += p / (self.sigma**k)
        return s * rising_factorial(self.theta + 1, n - 1) / rising_factorial(self.theta + 1, n + m)
