#= Module for plotting Pearson correlation coefficients ρ from timeseries data

!note: some functionality is limited (such as no diff. markers), as the code was originally 
       written using a specific personalized theme that is not included in this repository
=#
#/ Start module
module PearsonPlotter

#/ Packages
using CairoMakie
using LaTeXStrings

# using Distributions
using Statistics
# using Optim
# using SpecialFunctions

using JLD2
using CSV, DataFrames, DataFramesMeta
using FHist

#~ Specify paths
#!note: if these do not exist, create them
const RDATAPATH = "../data/rdata/"
const CSVDATAPATH = "../data/csv/"
const JLDATAPATH = "../data/jld/"

include("../juliacode/compute-moments.jl")
using .Moments

#################
### FUNCTIONS ###
"""
		Plot probability densities of fluctuations in rescaled log abundances
    Aggregates all environments into a single plot
"""
function plot_pearson(;
    prefix::String = "longitudinal/",
    envstatsfname::String = CSVDATAPATH * prefix * "environmentstats.csv",
    jlddir::String = JLDATAPATH * prefix,
    nbins = 51,
    savefig=false,
    figname=nothing
)
    #/ Create figure
    width = .9*246
    fig = Figure(;
        size=(width,width/1.67), figure_padding=(2,8,2,5), backgroundcolor=:transparent
    )
    ax = Axis(
        fig[1,1],
        limits=(-1.02,1.02,1e-4,1e1),
        xlabel=L"\textrm{Pearson}\;r", ylabel=L"\textrm{pdf}",
        xlabelsize=12, ylabelsize=12,
        yscale=log10,
        yminorticksvisible=false,
        xticklabelsize=9, yticklabelsize=9,
        # yticks=[1e-4,1e-2,1e0,1e1]
    )
    #/ Specify colors 
    colors = CairoMakie.to_colormap(:tab10)
    markers = [
        :circle, :utriangle, :cross, :rect, :diamond,
        :dtriangle, :pentagon, :xcross, :hexagon, :rtriangle
    ]

    #/ Load environment names
    edb = CSV.read(envstatsfname, DataFrame, delim=", ")
    #/ Filter envnames to include only those for which a histogram exists
    #~!note: these should total 9 distinct environments
    edb = filter(row -> isfile(jlddir*"pearsoncorrelation_$(row.environmentname).jld2"), edb)
    
    #/ Scatter plot the histograms of rescaled log frequencies
    for (i, envname) in enumerate(edb.environmentname)
        filename = JLDATAPATH * prefix * "pearsoncorrelation_$(envname).jld2"
        pearsondb = JLD2.load(filename)["pearson"]
        ρ = pearsondb.pearson
        fh = FHist.Hist1D(ρ; binedges=range(-1.,+1.,nbins), counttype=Int) |> normalize

        #@ WORK IN PROGRESS
        #~ Fit a distribution
        # function negloglikelihood(θ, data)
        #     ξ, ω, α = θ
        #     d = SkewNormal(ξ, ω^2, α)
        #     return -sum(logpdf.(d, data))
        # end
        # result = optimize(θ -> negloglikelihood(θ, ρ), [1.0, 1.0, 0.0])
        # θhat = Optim.minimizer(result)
        # ξ, ω, α = θhat
        # _fit = SkewNormal(ξ, ω^2, α)
        
        #/ Plot fit
        # fitplot = collect(range(-1., 1., 256))
        # lines!(
        #     ax, fitplot, Distributions.pdf.(_fit, fitplot),
        #     color=colors[i], linewidth=.8
        # )
        
        #~ Compute x-values at which to plot
        xplot = (fh.binedges[begin][2:end] + fh.binedges[begin][1:end-1]) ./ 2
        
        #/ Plot histograms
        vlines!(ax, [0.], color=:gray, linestyle=(:dash,:dense), linewidth=.6)
        scatter!(
            ax, xplot, fh.bincounts, markersize=3.5, strokewidth=.4,
            marker=markers[i], color=colors[i], label=envname
        )
    end
    
    #/ Add legend(s)
    Legend(
        fig[1,2], ax, labelsize=8, rowgap=0, patchsize=(2,2), framevisible=false
    )
    colgap!(fig.layout, 1)
    resize_to_layout!(fig)
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module MADPlotter
#/ End module
