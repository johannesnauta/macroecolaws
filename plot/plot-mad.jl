#= Module for plotting mean abundance distributions (MADs) from data

!note: some functionality is limited (such as no diff. markers), as the code was originally 
       written using a specific personalized theme that is not included in this repository
=#
#/ Start module
module MADPlotter

#/ Packages
using CairoMakie
using LaTeXStrings

using Distributions
using Statistics
using SpecialFunctions

using JLD2
using CSV, DataFrames, DataFramesMeta
using FHist

#~ Specify paths
#!note: if these do not exist, create them 
const RDATAPATH = "../data/rdata/"
const CSVDATAPATH = "../data/csv/"
const JLDATAPATH = "../data/jld/"

#################
### FUNCTIONS ###
"""
		Plot probability densities of rescaled log abundances
    Aggregates all environments into a single plot
"""
function plot_mad(;
    envstatsfname::String = CSVDATAPATH * "environmentstats.csv",
    histdir::String = JLDATAPATH,
    rescale=true,
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
        limits=(-6,6,1e-4,1e1),
        xlabel=L"\textrm{rescaled\;log\;abundances}", ylabel=L"\textrm{pdf}",
        xlabelsize=12, ylabelsize=12,
        yscale=log10, yminorticksvisible=false,
        xticklabelsize=9, yticklabelsize=9
    )    

    #/ Load environment names
    edb = CSV.read(envstatsfname, DataFrame, delim=", ")
    #/ Filter envnames to include only those for which a histogram exists
    #~!note: these should total 9 distinct environments
    edb = filter(row -> isfile(histdir*"fhist_$(row.environmentname).jld2"), edb)

    for (i, envname) in enumerate(edb.environmentname)
        #/ Load histogram and normalize
        fh = JLD2.load(histdir*"fhist_$(envname).jld2")["histogram"] |> normalize
        #~ Compute x-values at which to plot
        xplot = (fh.binedges[begin][2:end] + fh.binedges[begin][1:end-1]) ./ 2
        xpdf = range(-5,5,250)
        #/ Plot
        if !rescale
            #/ Just plot directly
            scatter!(ax, xplot, fh.bincounts, markersize=2, strokewidth=.5, label=envname)
        else
            #/ Rescale and then plot
            μ = edb[i,:mu]
            σ = edb[i,:sigma]
            c = edb[i,:cutoff]

            #/ Strange Grilli rescaling
            # (i == 1) && (
            #     lines!(ax, xpdf,10.0.^(-xpdf.^2), color=:black,linewidth=0.5,linestyle=:dash)
            # )
            # xscaled = @. (xplot - μ) / σ / sqrt(2)
            # pdfscaled = @. 10^(log(fh.bincounts)-log(0.5*erfc(μ-c)/sqrt(2*σ^2) + log(2π)/2))

            #/ Normal rescaling
            if i == 1
                n = Normal(0,1)
                lines!(ax, xpdf, xpdf -> pdf(n,xpdf), color=:black, linewidth=0.8)
            end
            xscaled = @. (xplot - μ) / σ
            pnormal = Normal(μ,σ)
            Z = 1 - cdf(pnormal, c)
            pdfscaled = fh.bincounts * σ * Z
            #~ Plot
            scatter!(ax, xscaled, pdfscaled, markersize=4, strokewidth=.5, label=envname)
        end
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

"""
    Plot MADs for each environment separately
    Does not aggregate, intended to display the differences between the environments, which
    is lost when one rescales the distributions
"""
function plot_mads(;
    envstatsfname::String = CSVDATAPATH * "environmentstats.csv",
    histdir::String = JLDATAPATH,
    rescale = false,
    savefig=false, figname=nothing
)
    #/ Create figure
    width = 1.4 * 246
    fig = Figure(;
        size=(width,width/1.5), figure_padding=(2,8,2,2), backgroundcolor=:transparent
    )
    #/ Specify colors 
    colors = CairoMakie.to_colormap(:tab10)
    markers = [
        :circle, :utriangle, :cross, :rect, :diamond,
        :dtriangle, :pentagon, :xcross, :hexagon, :rtriangle
    ]

    #/ Load environment names
    db = CSV.read(envstatsfname, DataFrame, delim=", ")
    #/ Filter envnames to include only those for which a histogram exists
    #~!note: these should total 9 distinct environments
    db = filter(row -> isfile(histdir*"fhist_$(row.environmentname).jld2"), db)

    ax = [Axis(
        fig[i,j],
        title=L"\textrm{%$(db.environmentname[3*(i-1)+j])}", titlesize=9, titlegap=1,
        xlabel=L"\textrm{log\;abundances}", xlabelvisible=(i==3),
        ylabel=L"\textrm{pdf}", ylabelvisible=(j==1),
        xlabelsize=10, ylabelsize=10,
        yscale=log10, yminorticksvisible=false,
        xticklabelsvisible=(i==3), yticklabelsvisible=(j==1),
        limits=(-25,0,1e-3,1e0),
        xticklabelsize=7, yticklabelsize=7,
    ) for i in 1:3, j in 1:3]
    
    # #~ Plot
    xpdf = range(-20,5,250)
    for (i, envname) in enumerate(db.environmentname)
        #/ Load histogram
        fh = JLD2.load(histdir*"fhist_$(envname).jld2")["histogram"] |> normalize
        
        #/ Compute grid indices
        gi = (i-1) ÷ 3 + 1
        gj = mod1(i,3)
        
        c = only(db[db.environmentname .== envname, :cutoff])
        xplot = (fh.binedges[begin][2:end] + fh.binedges[begin][1:end-1]) ./ 2

        if !rescale                      
            μ = db[i,:mu]
            σ = db[i,:sigma]
            c = db[i,:cutoff]
            n = truncated(Normal(μ, σ), lower=c)
            lines!(ax[gi,gj], xpdf, xpdf -> pdf(n,xpdf), color=:black, linewidth=0.5)
            
            vlines!(ax[gi,gj], [c], linestyle=(:dash,:dense), color=:black, linewidth=.5)
            scatter!(
                ax[gi,gj], xplot, fh.bincounts, markersize=4,
                strokewidth=.5, strokecolor=colors[i], marker=markers[i], color=:white,
                label=envname
            )
        else
            nothing
            #@TODO Implement rescaling
        end
    end
    #/ Add legend
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 10)
    resize_to_layout!(fig)
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module MADPlotter
#/ End module
