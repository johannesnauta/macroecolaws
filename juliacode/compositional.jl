#= Simple module that implements function to sample from gamma distributions and looks at
   other values that can be derived from them, like the distribution of the means, etc.
=#
#/ Start module
module CompSampler

using CairoMakie

using Distributions
using Random
using StatsBase

function sample_afd(
    S::Int;
    nsamples=512,
    seed::Int=42,
    rng=Random.Xoshiro(seed*S),
    dist=Distributions.Gamma
)
    if insupport(dist, -1)
        #~ truncate
        distf(α,θ) = truncated(dist(α,θ), lower=0.)
    else
        distf = dist
    end
    #/ Set shape- and scale-parameter(s)
    #~ note; for other distributions, they may not represent shape/scale
    α = rand(rng, S) .+ 5
    θ = 1.0 .+ 0.1 .* randn(rng, S)
    #/ Sample matrix and normalize so that ∑ₙxₙ=1
    x = reduce(vcat, [rand(rng, distf(α[i], θ[i]), nsamples) for i in eachindex(α)]')
    x = log.(x ./ sum(x, dims=1)) 
    #~ standardize
    μx = mean(x, dims=2)
    σx = std(x, dims=2)
    xrescaled = (x .- μx) ./ σx
    return x, xrescaled
end

function plot_afd(
    x::Matrix{Float64};
    nbins::Int = 11
)
    width = .8 * 246
    fig = Figure(; size=(width,width/1.33))
    ax = Axis(
        fig[1,1],
        xlabel="x", ylabel="pdf", yscale=log10,
        xlabelsize=10, ylabelsize=10, xticklabelsize=7, yticklabelsize=7
    )

    xplot = vec(x)
    xmin, xmax = extrema(xplot)
    bins = range(xmin, xmax, nbins)
    h = hist!(
        ax, xplot, bins=bins, normalization=:pdf,
        strokewidth=1, strokecolor=:black, color=:rebeccapurple
    )
    currentylims = ax.yaxis.attributes.limits.val
    #~ fit gamma distribution
    ex = exp.(xplot)
    shape = mean(ex).^2 / var(ex)
    scale = var(ex) / mean(ex)
    gammamoments = Distributions.Gamma(shape, scale)
    xfit = range(bins[begin], bins[end], 128)
    lines!(ax, xfit, exp.(xfit) .* Distributions.pdf.(gammamoments, exp.(xfit)), color=:red)
    # gammafit = Distributions.fit_mle(LogitNormal, ex)
    # lines!(ax, xfit, exp.(xfit) .* Distributions.pdf.(gammafit, exp.(xfit)), color=:black)
    return fig
end

end # module GSampler
#/ End module
