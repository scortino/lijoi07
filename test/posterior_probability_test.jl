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
x = 1:m

# Plots customization
plot_1 = plot(
    dpi = 300,
    size = (900, 600),
    left_margin = 7.5mm,
    xlabel = "\$k\$",
    ylabel = "pr\$(K_{50}^{(50)} = k | K_{50} = 5)\$",
    xlims = (0, 51),
    ylims = (0, 0.180),
    xticks = 0:5:50,
)
plot_2 = deepcopy(plot_1)
plot!(plot_2, ylabel = "pr\$(K_{50}^{(50)} = k | K_{50} = 25)\$")
plot_3 = deepcopy(plot_1)
plot!(plot_3, ylabel = "pr\$(K_{50}^{(50)} = k | K_{50} = 45)\$")

# Example 1: Dirichlet process
θ_dp = 19.233
p = DirichletProcess(θ_dp)
probs_dp = [posterior_probability(p, big(m), k, n) for k in x] # big necessary to compute factorial(21)
for plot_ in [plot_1, plot_2, plot_3]
    plot!(
        plot_,
        x,
        probs_dp,
        markershape = :xcross,
        markersize = 2,
        label = "DP \$(\\theta = 19.233)\$",
    )
end
@test sum(probs_dp) ≈ 1 atol = 1e-5

# Merge subplots
plot(plot_1, plot_2, plot_3, layout = (3, 1))

# Saving final plot
if !isfile("../img/posterior_probability.png")
    savefig("../img/posterior_probability.png")
end
