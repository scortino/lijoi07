from scipy.special import comb

from .utils import rising_factorial, generalized_factorial, noncentral_generalized_factorial, stirling

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

    def posterior_prob_K(self, m, k, n, j):
        s = 0
        for l in range(k, m + 1):
            s += comb(m, l) * stirling(l, k) * rising_factorial(n, m - l)
        return s * self.theta**k * rising_factorial(self.theta, n) / rising_factorial(self.theta, n + m)
    
    def Dhat_m(self, m, n, j):
        s1 = 0
        for k in range(m + 1):
            s2 = 0 
            for l in range(k, m + 1):
                s2 += comb(m, l) * stirling(l, k) * rising_factorial(n, m - l)
            s1 += self.theta**k * s2
        return s1 * self.theta / rising_factorial(self.theta + n, m + 1)
