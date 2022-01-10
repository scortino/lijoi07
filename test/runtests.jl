using Test

@testset "Function argument validation" begin
    include("validation.jl")
end

@testset "Compute prior probabilities for K_50" begin
    include("prior_probability_test.jl")
end

@testset "Compute posterior probabilities for K_50^(50)" begin
    include("posterior_probability_test.jl")
end
