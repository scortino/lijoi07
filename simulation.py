import matplotlib.pyplot as plt

from src.priors import DirichletProcess

### Section 4.1 A simple numerical example ###
n = 50
j = 25
m = 50

# Example 1: The Dirichlet Process
theta_dp = 19.233
dp = DirichletProcess(theta_dp)


# Plots
# Prior probabilities
prior_probs_dp = [dp.prior_prob_K(n, k) for k in range(1, n + 1)]

x = range(1, n + 1)
plt.figure()
plt.plot(x, prior_probs_dp, 'x-g', markersize=3, linewidth=0.5, label=r'DP$(\theta = 19.233)$')
plt.suptitle(r'Prior probabilities for $K_{50}$')
plt.xlabel(r'$k$')
plt.ylabel(r'pr$(K_{50} = k)$')
plt.axis([0, 51, 0, 0.125])
plt.xticks(range(0, 51, 5))
plt.legend()
plt.savefig('./img/prior_probabilities.png')
plt.close()
