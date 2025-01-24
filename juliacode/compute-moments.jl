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
"""Compute moments of a lower truncated lognormal distribution using MLE"""
function getestimates(logdata; c::Float64 = -100.0)
    #/ Compute empirical moments
    m1 = Statistics.mean(logdata)
    m2 = Statistics.mean(logdata.^2)

    #/ Define function that needs to be optimized
    function mle(dμ, μ, p)
        m1, m2, c = p
        σ = @. NaNMath.sqrt(-c*m1 + m2 + μ*(c - m1))
        x = @. (c - μ) / sqrt(2*σ^2)
        dμ .= @. (μ - m1)*erfc(x) - sqrt(2*σ^2 / π) * exp(-x^2)
        return nothing
    end

    nlprob = NonlinearProblem(mle, [m1], (m1, m2, c))
    #@TODO; specify a method that appears the most stable
    nlsol = solve(nlprob, TrustRegion(), reltol=1e-4, abstol=1e-4)
    μest = only(nlsol.u)
    σest = sqrt(-c*m1 + m2 + μest*(c - m1))
    return μest, σest
end


"""Compute moments of a truncated lognormal distribution"""
function fittrunclognormal(samples; uguess = [0.0, 1.0], lower=-Inf, upper=Inf)
    #/ 
    function truncnormlikelihood(p, data)
        μ, σ = p
        p = truncated(Normal(μ,sqrt(σ^2)), lower=lower, upper=upper)
        return -sum(logpdf.(p, data))
    end

    result = Optim.optimize(x -> truncnormlikelihood(x, samples), uguess, Optim.NelderMead())
    return Optim.minimizer(result)
end

end # module Moments
#/ End module
