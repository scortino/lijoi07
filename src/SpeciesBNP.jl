module SpeciesBNP

module Utils
include("utils.jl")
end

export 
    DirichletProcess,
    PoissonDirichletProcess,
    NormalizedIGProcess,
    prior_probability,
    posterior_probability

include("priors.jl")

end
