using Plots, Plots.PlotMeasures
using SpeciesBNP
using Test

# Note: when integers/floats casted to BigInt/BigFloat interact with others in operations,
# the latter are automatically promoted to BigInt/BigFloat
# https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic

### Section 4.1 A simple numerical example ###
n = 50
m = 50
js = [5, 25, 45]
x = 0:m # there can be 0 new species observed since we have already observed some

# Plots customization
plot_1 = plot(
    dpi = 300,
    size = (900, 600),
    left_margin = 7.5mm,
    xlabel = "\$k\$",
    ylabel = "pr\$(K_{50}^{(50)} = k | K_{50} = 5)\$",
    xlims = (-1, 51),
    ylims = (0, 0.175),
    xticks = 0:5:50,
    yticks = 0.025:0.025:0.175,
    guidefontsize = 9,
)
plot_2 = deepcopy(plot_1)
plot!(plot_2, ylabel = "pr\$(K_{50}^{(50)} = k | K_{50} = 25)\$")
plot_3 = deepcopy(plot_1)
plot!(plot_3, ylabel = "pr\$(K_{50}^{(50)} = k | K_{50} = 45)\$")
plots = [plot_1, plot_2, plot_3]

# Example 1: Dirichlet process
θ_dp = 19.233
p = DirichletProcess(θ_dp)
probs_dp = [posterior_probability(p, big(m), k, n) for k in x] # big necessary to compute factorial(21)
for plot_ in plots
    plot!(
        plot_,
        x,
        probs_dp,
        markershape = :xcross,
        markersize = 2,
        label = "DP \$(\\theta = 19.233)\$",
    )
end
@test sum(probs_dp) ≈ 1 atol = 1e-8

# Example 2: Two-parameter Poisson-Dirichlet process
σ_pd_1, θ_pd_1 = (0.25, 12.216)
p_1 = PoissonDirichletProcess(σ_pd_1, θ_pd_1)

σ_pd_2, θ_pd_2 = (0.75, 0.698)
p_2 = PoissonDirichletProcess(σ_pd_2, θ_pd_2)

for (i, j) in enumerate(js)
    local probs_pd_1 = [posterior_probability(p_1, big(m), big(k), n, j) for k in x] # big necessary to compute factorial(21)
    plot!(
        plots[i],
        x,
        probs_pd_1,
        markershape = :circle,
        markersize = 2,
        label = "PD \$(\\sigma = 0.25, \\theta = 12.216)\$",
    )
    @test sum(probs_pd_1) ≈ 1 atol = 1e-8

    local probs_pd_2 = [posterior_probability(p_2, big(m), big(k), n, j) for k in x] # big necessary to compute factorial(21)
    plot!(
        plots[i],
        x,
        probs_pd_2,
        markershape = :circle,
        markersize = 2,
        linestyle = :dash,
        label = "PD \$(\\sigma = 0.75, \\theta = 0.698)\$",
    )
    @test sum(probs_pd_2) ≈ 1 atol = 1e-4 # j = 45 most imprecise
end

# Example 3: normalized-inverse Gaussian process
x = 1:m # Posterior prob for N-IG defined only when k >= 1
θ_nig = 11.074
p = NormalizedIGProcess(big(θ_nig)) # big necessary to avoid approximation error

for (i, j) in enumerate(js)
    local probs_nig = [posterior_probability(p, big(m), k, big(n), j) for k in x] # big necessary to compute binomial(99, 18)
    plot!(
        plots[i],
        x,
        probs_nig,
        markershape = :utriangle,
        markersize = 2,
        linestyle = :dash,
        label = "N-IG \$(\\theta = 11.074)\$",
    )
    @test sum(probs_nig) ≈ 1 atol = 1e-4 # j = 5 most imprecise
end

# Merging subplots
plot(
    plot_1,
    plot_2,
    plot_3,
    layout = (3, 1),
    plot_title = "Posterior probability distributions for \$K_{50}^{(50)}\$",
)

# Saving final plot
if !isfile("../img/posterior_probability.png")
    savefig("../img/posterior_probability.png")
end
