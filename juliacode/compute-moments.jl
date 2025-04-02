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
"""Compute moments of a truncated lognormal distribution"""
function fittrunclognormal(samples; uguess = [-10.0, 1.0], lower=-Inf, upper=Inf)
    #/ Fit a truncated lognormal distribution
    function negtruncnormlikelihood(p, data)
        μ, logσ = p
        truncnorm = truncated(Normal(μ, exp(logσ)), lower=lower, upper=upper)
        return -sum(logpdf.(truncnorm, data))
    end

    result = Optim.optimize(
        x -> negtruncnormlikelihood(x, samples), uguess, Optim.NelderMead()
    )
    (Optim.converged(result)) && (return Optim.minimizer(result))
    @info "Optimizer not converged, returning guesses"
    return uguess
end

end # module Moments
#/ End module
