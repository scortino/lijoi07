using SpecialFunctions: gamma

include("utils.jl")

# /sigma -> 0
struct DirichletProcess
    θ::AbstractFloat
end

function prior_probability(prior::DirichletProcess, n::Integer, k::Integer)
    return prior.θ^k * stirling(n, k) / rising_factorial(prior.θ, n)
end

struct PoissonDirichletProcess
    σ::AbstractFloat
    θ::AbstractFloat
end

function prior_probability(prior::PoissonDirichletProcess, n::Integer, k::Integer)
    p = 1
    for i = 1:(k-1)
        p *= (prior.θ + i * prior.σ)
    end
    return p * generalized_factorial(n, k, prior.σ) / (prior.σ^k) /
           rising_factorial(prior.θ + 1, n - 1)
end

# /sigma = 1/2
struct NormalizedIGProcess
    θ::AbstractFloat
end

# When gamma takes two parameters, it evaluates to incomplete gamma function
function prior_probability(prior::NormalizedIGProcess, n::Integer, k::Integer)
    s = 0
    for i = 0:(n-1)
        s +=
            binomial(n - 1, i) * (-(prior.θ)^2)^(-i) * gamma(k + 2 + 2 * i - 2 * n, prior.θ)
    end
    return s * binomial(2 * n - k - 1, n - 1) * exp(prior.θ) * (-(prior.θ)^2)^(n - 1) /
           (2^(2 * n - k - 1) * gamma(k))
end
