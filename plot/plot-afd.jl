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

using CurveFit, NonlinearSolve

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
    envstatsfname::String = CSVDATAPATH * "environmentstats.csv",
    histdir::String = JLDATAPATH,
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
        limits=(-10,10,1e-4,1e0),
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
    edb = filter(row -> isfile(histdir*"afdfhist_$(row.environmentname).jld2"), edb)

    xfit = Float64[]
    yfit = Float64[]

    for (i, envname) in enumerate(edb.environmentname)
        #/ Load histogram and normalize
        fh = JLD2.load(histdir*"afdfhist_$(envname).jld2")["histogram"] |> normalize
        #~ Compute x-values at which to plot
        xplot = (fh.binedges[begin][2:end] + fh.binedges[begin][1:end-1]) ./ 2

        append!(xfit, exp.(xplot))
        append!(yfit, fh.bincounts)
        
        #/ Plot
        scatter!(
            ax, xplot, fh.bincounts, markersize=4, strokewidth=.5,
            marker=markers[i], color=colors[i], label=envname
        )
    end
    
    #/ Add legend
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
