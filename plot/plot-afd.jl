#= Module for plotting abundance fluctuation distributions (AFDs) from data

!note: some functionality is limited (such as no diff. markers), as the code was originally 
       written using a specific personalized theme that is not included in this repository
=#
#/ Start module
module AFDPlotter

#/ Packages
using CairoMakie
using LaTeXStrings

using Distributions
using Statistics
using SpecialFunctions

using JLD2
using CSV, DataFrames, DataFramesMeta
using FHist

using NonlinearSolve

#~ Specify paths
#!note: if these do not exist, create them
const RDATAPATH = "../data/rdata/"
const CSVDATAPATH = "../data/csv/"
const JLDATAPATH = "../data/jld/"

#################
### FUNCTIONS ###
"""
		Plot probability densities of fluctuations in rescaled log abundances
    Aggregates all environments into a single plot
"""
function plot_afd(;
    prefix::String = "longitudinal/",
    envstatsfname::String = CSVDATAPATH * prefix * "environmentstats.csv",
    jlddir::String = JLDATAPATH * prefix,
    compute_moments = true,
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
        limits=(-10,6,1e-4,1e0),
        xlabel=L"\textrm{rescaled\;log\;abundances}", ylabel=L"\textrm{pdf}",
        xlabelsize=12, ylabelsize=12,
        yscale=log10, yminorticksvisible=false,
        xticklabelsize=9, yticklabelsize=9
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
    edb = filter(row -> isfile(jlddir*"sampleafdfhist_$(row.environmentname).jld2"), edb)
    
    #/ Compute and save moments if they do not exist, otherwise load them
    if compute_moments
        #~ Collect frequencies of all environments
        freqs = Float64[]
        for (i, envname) in enumerate(edb.environmentname)
            filename = CSVDATAPATH * prefix * "rescaledlogfrequencydata_$(envname).csv"
            freqdb = CSV.read(filename, DataFrame, delim=", ")
            @info "hm" freqdb.sample_log_frequency
            append!(freqs, exp.(freqdb[!,:sample_log_frequency]))
        end
        @info "here?" freqs
        #~ Fit gamma distribution
        gammafit = Distributions.fit_mle(Gamma, freqs)
        lognormfit = Distributions.fit_mle(LogNormal, freqs)
        α, θ = params(gammafit)
        μ, σ = params(lognormfit)
        jldsave(jlddir * "gammaparams.jld2"; α=α, θ=θ)
        jldsave(jlddir * "lognormalparams.jld2"; μ=μ, σ=σ)
    else
        gammaparams = JLD2.load(jlddir * "gammaparams.jld2")
        α, θ = gammaparams["α"], gammaparams["θ"]
        lognormparams = JLD2.load(jlddir * "lognormalparams.jld2")
        μ, σ = lognormparams["μ"], lognormparams["σ"]
    end
    
    #/ Scatter plot the histograms of rescaled log frequencies
    for (i, envname) in enumerate(edb.environmentname)
        filename = CSVDATAPATH * prefix * "rescaledlogfrequencydata_$(envname).csv"
        freqdb = CSV.read(filename, DataFrame, delim=", ")
        gammafit = Distributions.fit_mle(Gamma, exp.(freqdb.log_frequency))
        #/ Load histogram and normalize
        fh = JLD2.load(jlddir*"sampleafdfhist_$(envname).jld2")["histogram"] |> normalize
        #~ Compute x-values at which to plot
        xplot = (fh.binedges[begin][2:end] + fh.binedges[begin][1:end-1]) ./ 2
        #~ Rescale to standard gamma
        # xplot = xplot .- (gammafit.α * gammafit.θ)
        # yplot = fh.bincounts .* sqrt(gammafit.α * gammafit.θ^2)        
        
        #/ Plot
        scatter!(
            ax, xplot, fh.bincounts, markersize=3.5, strokewidth=.4,
            marker=markers[i], color=colors[i], label=envname
        )
    end  

    #/ Plot fitted gamma distribution
    xfits = exp.(range(-9, 5, 256))
    ygamma = xfits .* Distributions.pdf.(Gamma(α,θ), xfits)
    ylognormal = xfits .* Distributions.pdf(LogNormal(μ,σ), xfits)
    gammaline = lines!(ax, log.(xfits), ygamma, color=:black, linewidth=1.)
    lognormalline = lines!(
        ax, log.(xfits), ylognormal, color=:black, linewidth=.8, linestyle=(:dash,:dense)
    ) 
    
    #/ Add legend(s)
    axislegend(
        ax, [gammaline, lognormalline], [L"\textrm{gamma}", L"\textrm{lognormal}"],
        position=:lt, labelsize=8, nbanks=1, patchlabelgap=1.2,
        patchsize=(6,1), padding=0, margin=(2,0,0,2), framevisible=false
    )
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
