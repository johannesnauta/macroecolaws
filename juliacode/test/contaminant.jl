#/ Start module
module CTest

using Distributions, Optim, StatsBase, Plots, Random

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
    return Optim.minimizer(mle)
end

function plot_stuff()
    Random.seed!(42)
    n = 512
    ε = 0.97
    nε = round(Int, n * ε)
    samples = vcat(rand(Gamma(3, 2), nε), rand(Normal(25., 1.), n-nε))
    fitresult = fit_mixture(samples, prior=Distributions.Gamma)
    logshape, logscale, sigε, μ, logσ = fitresult
    ε = 1 / (1 + exp(-sigε))

    plt = histogram(samples, nbins=101, normalize=:pdf, legend=true)
    Γ = Gamma(exp(logshape), exp(logscale))
    
    G = Normal(μ, exp(logσ))
    xplot = range(1e-3, 1.25*maximum(samples), 256)
    Γpdf = pdf.(Γ, xplot)
    Gpdf = pdf.(G, xplot)
    plot!(xplot, Γpdf, color=:rebeccapurple, linestyle=:dash, label="gamma")
    plot!(xplot, Gpdf, color=:firebrick2, linestyle=:dash, label="normal")
    P = @. ε*Γpdf + (1-ε)*Gpdf
    plot!(plt, xplot, P, color=:black, label="mixture")
    return plt
end

# # Example Usage
# data = vcat(rand(Gamma(3, 2), 1000), rand(Uniform(0, 50), 50))  # Mostly Gamma, some outliers
# α_hat, θ_hat, ε_hat = fit_gamma_mixture(data)

# println("Fitted Gamma Shape: ", α_hat)
# println("Fitted Gamma Scale: ", θ_hat)
# println("Estimated Contamination Fraction: ", ε_hat)


end # module CTest
#/ End module
