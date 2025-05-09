#= Module for plotting Pearson correlation coefficients ρ from timeseries data

!note: some functionality is limited (such as no diff. markers), as the code was originally 
       written using a specific personalized theme that is not included in this repository
=#
#/ Start module
module PearsonPlotter

#/ Packages
using CairoMakie
using LaTeXStrings

using StatsBase
using Distributions
using Optim
using SpecialFunctions

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
    fitdistribution = false,
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
        limits=(-7.02,10.02,1e-5,1e1),
        xlabel=L"\textrm{rescaled Fisher}\;z", ylabel=L"\textrm{pdf}",
        xlabelsize=11, ylabelsize=11,
        yscale=log10,
        yminorticksvisible=false,
        xticklabelsize=9, yticklabelsize=9,
        yticks = LogTicks([-5,-3,-1,1]),
        xticks = [-5,0,5,10],
        # xgridvisible=false, ygridvisible=false
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

    rescaledz = []
    _mean = 0.
    _mode = 0.
    
    #/ Scatter plot the histograms of rescaled log frequencies
    for (i, envname) in enumerate(edb.environmentname)
        filename = JLDATAPATH * prefix * "pearsoncorrelation_$(envname).jld2"
        pearsondb = JLD2.load(filename)["pearson"]
        ρ = pearsondb.pearson
        # fh = FHist.Hist1D(ρ; binedges=range(-1.,+1.,nbins), counttype=Int) |> normalize
        zv = @. log((1+ρ)/(1-ρ)) / 2
        z = (zv .- mean(zv)) ./ std(zv)
        append!(rescaledz, z)
        bins = range(minimum(z),maximum(z),nbins)
        fh = FHist.Hist1D(z; binedges=bins, counttype=Int, overflow=true) |> normalize
        
        #@ WORK IN PROGRESS
        if i == length(edb.environmentname)
            fitplot = collect(range(-10., 10., 256))
            if fitdistribution
                rescaledz = dropdims(reduce(hcat, rescaledz), dims=1)
                #~ Fit a skewnormal distribution onto all the rescaled datapoints
                function negloglikelihood(θ, data)
                    ξ, ω, α = θ
                    d = SkewNormal(ξ, ω^2, α)
                    return -sum(logpdf.(d, data))
                end
                result = optimize(θ -> negloglikelihood(θ, rescaledz), [1.0, 1.0, 0.0])
                θhat = Optim.minimizer(result)
                ξ, ω, α = θhat
                _fit = SkewNormal(ξ, ω^2, α)
                _pdf = Distributions.pdf.(_fit, fitplot)
            else
                #~ Just plot the "standard" skewnormal distribution with mean 0, variance 1,
                #  and skewness of 1/2
                #~ note: when fitting it gives the same values, essentially, but as of writing I
                #        have no reason why it should be this skewed or something
                skewnorm = SkewNormal(-1,sqrt(2),2)
                _mode = StatsBase.mode(skewnorm)
                _mean = StatsBase.mean(skewnorm)
                _pdf = Distributions.pdf.(skewnorm, fitplot)
            end                
        
            #/ Plot fit
            lines!(ax, fitplot, _pdf, color=:black, linewidth=1.)
        end
        
        #~ Compute x-values at which to plot
        xplot = (fh.binedges[begin][2:end] + fh.binedges[begin][1:end-1]) ./ 2
        scatter!(
            ax, xplot, fh.bincounts, markersize=3.5, strokewidth=.4,
            marker=markers[i], color=colors[i], label=envname
        )
    end
            
    #/ Plot mode and mean    
    modeline = vlines!(ax, [_mode], color=:black, linestyle=(:solid,:dense), linewidth=.6)
    meanline = vlines!(ax, [_mean], color=:black, linestyle=(:dash,:dense), linewidth=.5)
    
    #/ Add legend(s)
    Legend(fig[1,2], ax, labelsize=8, rowgap=0, patchsize=(2,2), framevisible=false)
    axislegend(
        ax, [meanline, modeline], [L"\textrm{mean}", L"\textrm{mode}"],
        position=:rt, margin=(0,0,0,0),
        framevisible=false, labelsize=9, rowgap=0, patchsize=(8,1), 
    )
    colgap!(fig.layout, 1)
    resize_to_layout!(fig)
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module MADPlotter
#/ End module
