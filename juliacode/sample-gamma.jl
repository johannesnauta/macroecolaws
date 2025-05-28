#= Simple module that implements function to sample from gamma distributions and looks at
   other values that can be derived from them, like the distribution of the means, etc.
=#
#/ Start module
module GammaSampler

using Distributions
# using DataFrames, DataFramesMeta
using Random
using StatsBase
using NaNMath

function sample_afd(S::Int; nsamples=100, seed::Int=42, rng=Random.Xoshiro(seed*S))
    α = 1.0 .* ones(S)
    θpdf = LogNormal(1e-3, 1.)
    θ = rand(rng, θpdf, S)

    #~ Sample from the AFD
    x = zeros(S, nsamples)
    for i in 1:S
        x[i,:] = rand(rng, Gamma(α[i],θ[i]), nsamples)
    end
    return x
end

function rescaledpdf(x; cutoff=-Inf)
    #~ apply cutoff
    logx = log.(x)
    logx[x .< cutoff] .= NaN
    meanx = mapslices(x -> NaNMath.mean(x), logx; dims=2)
    stdx = mapslices(x -> NaNMath.std(x), logx; dims=2)
    logx = (logx .- meanx) ./ stdx
    return logx, meanx
end

function sample_mad(x; cutoff=-Inf)
    logx = rescaledpdf(x, cutoff=cutoff)
    meanx = mapslices(z -> NaNMath.mean(z), logx; dims=2)
    return logx, meanx
end

end # module GSampler
#/ End module
