using SpeciesBNP
using Test

# priors.py
p = DirichletProcess(1.0)
@test_throws DomainError prior_probability(p, 0, 0)
@test_throws DomainError prior_probability(p, 5, 0)
@test_throws DomainError prior_probability(p, 5, 6)
@test_throws DomainError posterior_probability(p, 0, 0, 10)
@test_throws DomainError posterior_probability(p, 5, -1, 10)
@test_throws DomainError posterior_probability(p, 5, 6, 10)
@test_throws DomainError posterior_probability(p, 5, 5, 0)

p = PoissonDirichletProcess(1.0, 1.0)
@test_throws DomainError prior_probability(p, 0, 0)
@test_throws DomainError prior_probability(p, 5, 0)
@test_throws DomainError prior_probability(p, 5, 6)
@test_throws DomainError posterior_probability(p, 0, 0, 10, 10)
@test_throws DomainError posterior_probability(p, 5, -1, 10, 10)
@test_throws DomainError posterior_probability(p, 5, 6, 10, 10)
@test_throws DomainError posterior_probability(p, 5, 5, 0, 10)
@test_throws DomainError posterior_probability(p, 5, 5, 10, 11)
@test_throws DomainError posterior_probability(p, 5, 5, 10, 0)

p = NormalizedIGProcess(1.0)
@test_throws DomainError prior_probability(p, 0, 0)
@test_throws DomainError prior_probability(p, 5, 0)
@test_throws DomainError prior_probability(p, 5, 6)
@test_throws DomainError posterior_probability(p, 0, 0, 10, 10)
@test_throws DomainError posterior_probability(p, 5, 0, 10, 10)
@test_throws DomainError posterior_probability(p, 5, 6, 10, 10)
@test_throws DomainError posterior_probability(p, 5, 5, 0, 10)
@test_throws DomainError posterior_probability(p, 5, 5, 10, 11)
@test_throws DomainError posterior_probability(p, 5, 5, 10, 0)

# utils.py
@test_throws DomainError SpeciesBNP.Utils.rising_factorial(1, -1)
@test_throws DomainError SpeciesBNP.Utils.generalized_factorial(-1, 0, 1.0)
@test_throws DomainError SpeciesBNP.Utils.generalized_factorial(0, -1, 1.0)
@test_throws DomainError SpeciesBNP.Utils.noncentral_generalized_factorial(-1, 0, 1.0, 1.0)
@test_throws DomainError SpeciesBNP.Utils.noncentral_generalized_factorial(0, -1, 1.0, 1.0)
@test_throws DomainError SpeciesBNP.Utils.stirling(-1, 0)
@test_throws DomainError SpeciesBNP.Utils.stirling(0, -1)
