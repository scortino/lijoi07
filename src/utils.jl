using LRUCache
using Memoize
using SpecialFunctions

# Rising factorial
function rising_factorial(a::Real, n::Integer)
    prod(a + j for j = 0:(n-1); init = 1)
end

# Central generalized factorial coefficient
function generalized_factorial(n::Integer, k::Integer, σ::AbstractFloat)
    if n == k == 0
        return 1
    end

    sum(
        (-1)^j * rising_factorial(-j * σ, n) / (factorial(j) * factorial(k - j)) for j = 0:k;
        init = 0,
    )
end

# Non-central generalized factorial coefficient
function noncentral_generalized_factorial(
    n::Integer,
    k::Integer,
    σ::AbstractFloat,
    γ::AbstractFloat,
)
    sum(
        binomial(n, j) * generalized_factorial(j, k, σ) * rising_factorial(-γ, n - j) for
        j = k:n;
        init = 0,
    )
end

# Unsigned stirling number of the first kind
# https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
# n>=0, k>=0
@memoize LRU{Tuple{Any,Any},Any}(maxsize = 128) function stirling(n::Integer, k::Integer)
    if n == k == 0
        return 1
    elseif (n == 0) || (k == 0)
        return 0
    elseif n == k
        return 1
    elseif k == 1
        return factorial(n - 1)
    elseif k == n - 1
        return binomial(n, 2)
    elseif k == n - 2
        return div((3 * n - 1) * binomial(n, 3), 4)
    elseif k == n - 3
        return binomial(n, 2) * binomial(n, 4)
    end

    (n - 1) * stirling(n - 1, k) + stirling(n - 1, k - 1)
end
