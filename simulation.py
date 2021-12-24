import matplotlib.pyplot as plt

from src.priors import DirichletProcess, PoissonDirichletProcess

### Section 4.1 A simple numerical example ###
n = 50
j = 25
m = 50

# Example 1: The Dirichlet Process
theta_dp = 19.233
dp = DirichletProcess(theta_dp)

# Example 2: The two-parameter Poisson-Dirichlet process
sigma_pd_1, theta_pd_1 = (0.25, 12.216)
sigma_pd_2, theta_pd_2 = (0.75, 0.698)
pd_1 = PoissonDirichletProcess(sigma_pd_1, theta_pd_1)
pd_2 = PoissonDirichletProcess(sigma_pd_2, theta_pd_2)


# Plots
# Prior probabilities
prior_probs_dp = [dp.prior_prob_K(n, k) for k in range(1, n + 1)]
prior_probs_pd_1 = [pd_1.prior_prob_K(n, k) for k in range(1, n + 1)]
prior_probs_pd_2 = [pd_2.prior_prob_K(n, k) for k in range(1, n + 1)]

x = range(1, n + 1)
plt.figure(figsize=(12, 8))
plt.plot(x, prior_probs_dp, 'x-g', markersize=3, linewidth=0.5, label=r'DP$(\theta = 19.233)$')
plt.plot(x, prior_probs_pd_1, '.-r', markersize=3, linewidth=0.5, label=r'PD$(\sigma = 0.25, \theta = 12.216)$')
plt.plot(x, prior_probs_pd_2, '.--b', markersize=3, linewidth=0.5, label=r'PD$(\sigma = 0.75, \theta = 0.698)$')
plt.suptitle(r'Prior probabilities for $K_{50}$')
plt.xlabel(r'$k$')
plt.ylabel(r'pr$(K_{50} = k)$')
plt.axis([0, 51, 0, 0.125])
plt.xticks(range(0, 51, 5))
plt.legend()
plt.tight_layout()
plt.savefig('./img/prior_probabilities.png')
plt.close()
