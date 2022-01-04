using Plots, Plots.PlotMeasures
using SpeciesBNP
using Test

# Note: when integers/floats casted to BigInt/BigFloat interact with others in operations,
# the latter are automatically promoted to BigInt/BigFloat
# https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic

### Section 4.1 A simple numerical example ###
n = 50
x = 1:n

# Example 1: Dirichlet process
θ_dp = 19.233
p = DirichletProcess(θ_dp)
probs_dp = [prior_probability(p, big(n), k) for k in x] # big necessary to compute factorial(49)
plot!(
    x,
    probs_dp,
    markershape = :xcross,
    markersize = 2,
    label = "DP \$(\\theta = 19.233)\$",
)
@test sum(probs_dp) ≈ 1 atol = 1e-8

# Example 2: Two-parameter Poisson-Dirichlet process
σ_pd_1, θ_pd_1 = (0.25, 12.216)
p = PoissonDirichletProcess(σ_pd_1, θ_pd_1)
probs_pd_1 = [prior_probability(p, n, big(k)) for k in x] # big necessary to compute factorial(21)
plot!(
    x,
    probs_pd_1,
    markershape = :circle,
    markersize = 2,
    label = "PD \$(\\sigma = 0.25, \\theta = 12.216)\$",
)
@test sum(probs_pd_1) ≈ 1 atol = 1e-8

σ_pd_2, θ_pd_2 = (0.75, 0.698)
p = PoissonDirichletProcess(σ_pd_2, θ_pd_2)
probs_pd_2 = [prior_probability(p, n, big(k)) for k in x] # big necessary to compute factorial(21)
plot!(
    x,
    probs_pd_2,
    markershape = :circle,
    markersize = 2,
    linestyle = :dash,
    label = "PD \$(\\sigma = 0.75, \\theta = 0.698)\$",
)
@test sum(probs_pd_2) ≈ 1 atol = 1e-8

# Example 3: normalized-inverse Gaussian process
θ_nig = 11.074
p = NormalizedIGProcess(big(θ_nig)) # big necessary to avoid approximation error
probs_nig = [prior_probability(p, n, big(k)) for k in x] # big necessary to compute binomial(98, 49)
plot!(
    x,
    probs_nig,
    markershape = :utriangle,
    markersize = 2,
    linestyle = :dash,
    label = "N-IG \$(\\theta = 11.074)\$",
)
@test sum(probs_nig) ≈ 1 atol = 1e-8

# Plot customization
plot!(
    dpi = 300,
    size = (900, 600),
    left_margin = 7.5mm,
    title = "Prior probabilities for \$K_{50}\$",
    xlabel = "\$k\$",
    ylabel = "pr\$(K_{50} = k)\$",
    xlims = (0, 51),
    ylims = (0, 0.125),
    xticks = 0:5:50,
)

# Saving final plot
if !isfile("../img/prior_probability.png")
    savefig("../img/prior_probability.png")
end
