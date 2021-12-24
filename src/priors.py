from utils import generalized_factorial, noncentral_generalized_factorial

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