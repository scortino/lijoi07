using SpecialFunctions: gamma

include("utils.jl")

# /sigma -> 0
struct DirichletProcess
    θ::AbstractFloat
end

function prior_probability(prior::DirichletProcess, n::Integer, k::Integer)
    return prior.θ^k * stirling(n, k) / rising_factorial(prior.θ, n)
end

# Notice it doesn't depend on j
function posterior_probability(prior::DirichletProcess, m::Integer, k::Integer, n::Integer)
    s = 0
    for l = k:m
        s += binomial(m, l) * stirling(l, k) * rising_factorial(n, m - l)
    end
    return s * prior.θ^k * rising_factorial(prior.θ, n) / rising_factorial(prior.θ, n + m)
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

function posterior_probability(
    prior::PoissonDirichletProcess,
    m::Integer,
    k::Integer,
    n::Integer,
    j::Integer,
)
    s = 1
    for i = j:(j+k-1)
        s *= (prior.θ + i * prior.σ)
    end
    return s / prior.σ^k *
           noncentral_generalized_factorial(m, k, prior.σ, -n + j * prior.σ) *
           rising_factorial(prior.θ + 1, n - 1) / rising_factorial(prior.θ + 1, n + m - 1)
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
