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
        return -sum(logpdf.(truncnorm, data[data.>lower]))
    end

    result = Optim.optimize(x -> neglikelihood(x, samples), uguess, Optim.NelderMead())
    (Optim.converged(result)) && (return Optim.minimizer(result))
    @info "Optimizer not converged, returning guesses"
    return uguess
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

end # module Moments
#/ End module
