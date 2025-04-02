#/ Start module
module Mixture

#/ Packages
using Distributions
using Optim

#~ Define struct to return
struct MLEStruct
    shape::Float64
    scale::Float64
    ε::Float64
    μ::Float64
    σ::Float64
end

# Define the mixture log-likelihood function
function negloglikelihood(params, data; prior = Distributions.Gamma)
    logα, logθ, sigε, μ, logσ = params  # Shape (α), scale (θ), and contamination fraction (ε)
    #/ Do the transformations
    #~ these transformations are done such that the parameters are in the correct domain
    #~ for ε we use a sigmoid as it's in (0,1)
    #~ for all other positive parameters, a simple log/exp transform
    ε = 1 / (1 + exp(-sigε))
    α = exp(logα)
    θ = exp(logθ)
    σ = exp(logσ)
    #/ Define the distributions, and compute their pdf
    p = prior(α, θ)
    N = Distributions.Normal(μ, σ)
    _priorpdf = pdf.(p, data)
    _normalpdf = pdf.(N, data)
    #/ Use the pdf ration p(x|ϑ) as the 'weight' for the data
    γ = @. ε * _priorpdf  / (ε * _priorpdf  + (1 - ε) * _normalpdf)
    loglikelihood = @. γ * log(ε * _priorpdf) + (1-γ)*log((1-ε) * _normalpdf)
    return -sum(loglikelihood)
end

# Fit the Gamma + Contaminant Model
function fit_mixture(samples; prior = Distributions.Gamma, guess = [1., 1., 0.5, 1., 1.])
    mle = optimize(x -> negloglikelihood(x, samples, prior=prior), guess, Optim.NelderMead())
    params = Optim.minimizer(mle)
    shape = exp(getindex(params, 1))
    scale = exp(getindex(params, 2))
    ε = 1 / (1 + exp(-getindex(params, 3)))
    μ = getindex(params, 4)
    σ = exp(getindex(params, 3))
    return MLEStruct(shape, scale, ε, μ, σ)
end

end # module Mixture
#/ End module
