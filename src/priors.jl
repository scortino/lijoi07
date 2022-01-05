using SpecialFunctions: gamma

include("utils.jl")

# /sigma -> 0
struct DirichletProcess
    θ::AbstractFloat
end

function prior_probability(prior::DirichletProcess, n::Integer, k::Integer)
    prior.θ^k * stirling(n, k) / rising_factorial(prior.θ, n)
end

# Notice it doesn't depend on j
function posterior_probability(prior::DirichletProcess, m::Integer, k::Integer, n::Integer)
    s = sum(
        binomial(m, l) * stirling(l, k) * rising_factorial(n, m - l) for l = k:m;
        init = 0,
    )
    prior.θ^k * rising_factorial(prior.θ, n) / rising_factorial(prior.θ, n + m) * s
end

struct PoissonDirichletProcess
    σ::AbstractFloat
    θ::AbstractFloat
end

function prior_probability(prior::PoissonDirichletProcess, n::Integer, k::Integer)
    p = prod(prior.θ + i * prior.σ for i = 1:(k-1); init = 1)
    p / (prior.σ^k * rising_factorial(prior.θ + 1, n - 1)) *
    generalized_factorial(n, k, prior.σ)
end

function posterior_probability(
    prior::PoissonDirichletProcess,
    m::Integer,
    k::Integer,
    n::Integer,
    j::Integer,
)
    p = prod(prior.θ + i * prior.σ for i = j:(j+k-1); init = 1)
    rising_factorial(prior.θ + 1, n - 1) / rising_factorial(prior.θ + 1, n + m - 1) * p /
    prior.σ^k * noncentral_generalized_factorial(m, k, prior.σ, -n + j * prior.σ)
end

# /sigma = 1/2
struct NormalizedIGProcess
    θ::AbstractFloat
end

# When gamma takes two parameters, it evaluates to incomplete gamma function
function prior_probability(prior::NormalizedIGProcess, n::Integer, k::Integer)
    s = sum(
        binomial(n - 1, i) * (-prior.θ^2)^(-i) * gamma(k + 2 + 2i - 2n, prior.θ) for
        i = 0:(n-1);
        init = 0,
    )
    binomial(2n - k - 1, n - 1) * (exp(prior.θ) * (-prior.θ^2)^(n - 1)) /
    (2^(2n - k - 1) * gamma(k)) * s
end

# Notice k >= 1
function posterior_probability(
    prior::NormalizedIGProcess,
    m::Integer,
    k::Integer,
    n::Integer,
    j::Integer,
)
    s1 = s2 = s3 = 0
    for i = 0:(n+m-1)
        s1 +=
            binomial(n + m - 1, i) *
            (-prior.θ^2)^(-i) *
            gamma(j + k + 2 + 2i - 2(m + n), prior.θ)
        if i < n
            s2 += binomial(n - 1, i) * (-prior.θ^2)^(-i) * gamma(j + 2 + 2i - 2n, prior.θ)
        end
        if k <= i <= m
            s3 +=
                binomial(m, i) * binomial(2i - k - 1, i - 1) * gamma(i) / 2^(2i) *
                rising_factorial(n - j / 2, m - i)
        end
    end
    ((-prior.θ^2)^m * 2^k) / (rising_factorial(n, m) * gamma(k)) * s1 / s2 * s3
end
