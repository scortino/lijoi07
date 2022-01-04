using Test

@testset "Compute prior probabilities for K_50" begin
    include("prior_probability_test.jl")
end

@testset "Compute posterior probabilities for K_50^(50)" begin
    include("posterior_probability_test.jl")
end
