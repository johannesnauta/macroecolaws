#/ Start module
module CTest

using Distributions, Optim, StatsBase, Plots, Random

# Define the mixture log-likelihood function
function gamma_mixture_loglik(params, data)
    logα, logθ, ε, μ, logσ = params  # Shape (α), scale (θ), and contamination fraction (ε)
    if ε < 0 || ε > 1
        ε = 1 / ε
    end
    Γ = Gamma(exp(logα), exp(logθ))
    G = Normal(μ, exp(logσ))
    Γpdf = pdf.(Γ, data)
    Gpdf = pdf.(G, data)

    γ = @. ε * Γpdf  / (ε * Γpdf  + (1 - ε) * Gpdf)
    Q = @. γ * log(ε * Γpdf) + (1-γ)*log((1-ε) * Gpdf)
    loglikelihood = sum(Q)
    return -loglikelihood  # Negative log-likelihood for minimization
end

# Fit the Gamma + Contaminant Model
function fit_gamma_mixture(data)
    MLE = fit_mle(Gamma, data)  # Initial estimates from standard MLE
    αguess, θguess = MLE.α, MLE.θ
    εguess = 0.5
    params = [log(αguess), log(θguess), εguess, maximum(data), log(std(data))]

    mixturefit = optimize(x -> gamma_mixture_loglik(x, data), params, Optim.NelderMead())
    return Optim.minimizer(mixturefit)
end

function plot_stuff()
    Random.seed!(42)
    samples = vcat(rand(Gamma(3, 2), 500), rand(Normal(100., 10.), 1))
    fitresult = fit_gamma_mixture(samples)
    logshape, logscale, ε, μ, logσ = fitresult

    plt = histogram(samples, nbins=101, normalize=:pdf, legend=true)
    Γ = Gamma(exp(logshape), exp(logscale))

    # U = Uniform(umin, umax)
    # U = Uniform(minimum(samples) + log(umin), maximum(samples) + log(umax))
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
