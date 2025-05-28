#= Module for computing moments a given distribution provided a dataset of frequencies =#
#/ Start module
module Moments

#/ Packages
using NaNMath
using NonlinearSolve
using SpecialFunctions
using Statistics

using Distributions
using Optim
using Random

#################
### FUNCTIONS ###
"""Compute parameters of a truncated lognormal distribution"""
function fittrunclognormal(samples; uguess = [-10.0, 1.0], lower=-Inf, upper=Inf)
    #/ Fit a truncated lognormal distribution
    function neglikelihood(p, data)
        μ, logσ = p
        truncnorm = truncated(Normal(μ, exp(logσ)), lower=lower, upper=upper)
        #~ make sure to check for data above the lower cutoff, otherwise the
        #  logpdf will have zeros leading to -Inf loglikelihoods
        return -sum(logpdf.(truncnorm, data[data .> lower]))
    end

    result = Optim.optimize(x -> neglikelihood(x, samples), uguess, Optim.NelderMead())
    (Optim.converged(result)) && (return Optim.minimizer(result))
    @info "Optimizer not converged, returning guesses"
    return uguess
end

"[wip] Nonlinearsolve for the MLE from Grilli"
function compute_MAD_params_nlsolve(m1, m2, c)
    function f!(dx, x, p)
        dx[1] = x[1] - m1 +
                sqrt(2/π) * x[2] * exp(-(c - x[1])^2 / (2 * x[2]^2)) /
                erfc((c - x[1]) / sqrt(2 * x[2]^2))
        dx[2] = x[2]^2 + m1*x[1] + c*m1 - x[1]*c - m2
    end

    prob = NonlinearProblem(f!, [-15.0, 2.0])
    sol = solve(prob, NewtonRaphson())

    return sol.u[1], sol.u[2]
end

"""Compute parameters of a generalized gamma distribution"""
function generalizedgamma(x; a=a, d=d, p=p)
    (x <= 0.0) && (return 0.0)
    return (p / (a^d)) * x^(d-1) * exp(-(x/a)^p) / Distributions.gamma(d/p)
end

function fitgeneralizedgamma(samples; uguess = [log(5.),log(2.),log(1.)])
	  function negloglikelihood(params, data)
	      loga, logd, logp = params
        ggamma = generalizedgamma.(data, a=exp(loga), d=exp(logd), p=exp(logp))
        return -sum(log.(ggamma))
    end

    result = Optim.optimize(x -> negloglikelihood(x, samples), uguess, Optim.NelderMead())
    (Optim.converged(result)) && (return Optim.minimizer(result))
    @info "Optimizer not converged, returning guesses"
    return uguess
end

function fit_mixture(data; uguess = [1e-8,1.,2.,1.,0.5])
    function negloglikelihood(p, data)
        α, θ, μ, logσ, w1 = p
        α = α^2
        θ = θ^2
        w = [w1^2, 1 - w1^2]
        model = MixtureModel([Gamma(α,θ), LogNormal(μ,exp(logσ))], w)
        return -sum(logpdf.(model, data))
    end

    result = Optim.optimize(x -> negloglikelihood(x, data), uguess, Optim.NelderMead())
    if Optim.converged(result)
        result = Optim.minimizer(result)
        α = result[1]^2
        θ = result[2]^2
        μ = result[3]
        σ = exp(result[4])
        w1 = result[5]^2
        w2 = 1 - w1
        return (; α=α, θ=θ, μ=μ, σ=σ, w=[w1,w2])
    end
    @info "Optimizer not converged, returning guesses"
    return uguess
end

end # module Moments
#/ End module
