#= Module for plotting mean abundance distributions (MADs) from data

!note: some functionality is limited (such as no diff. markers), as the code was originally 
       written using a specific personalized theme that is not included in this repository
=#
#/ Start module
module TaylorPlotter

#/ Packages
using CairoMakie
using LaTeXStrings

using CSV, JLD2
using CategoricalArrays
using DataFrames, DataFramesMeta
using Statistics

using CurveFit
using NonlinearSolve

#~ Specify paths
const RDATAPATH = "../data/rdata/"
const CSVDATAPATH = "../data/csv/"
const JLDATAPATH = "../data/jld/"

#################
### FUNCTIONS ###
"""
    Plot Taylor's law
    Aggregates all environments into a single plot
"""
function plot_taylor(;
    envstatsfname::String = CSVDATAPATH * "environmentstats.csv",
    nbins::Int = 15, raw=false, γ=2.0,
    freqdatadirname::String = dirname(envstatsfname) * "/",
    savefig=false, figname=nothing
)
    #/ Create figure
    width = .9*246
    fig = Figure(;
        size=(width,width/1.67), figure_padding=(2,8,2,5), backgroundcolor=:transparent
    )
    ax = Axis(
        fig[1,1],
        limits=(1e-9,1e0,1e-16,1e0),
        xlabel=L"\textrm{mean\;rel.\;abundance}",
        ylabel=L"\textrm{var.\;rel.\;abundance}",
        xlabelsize=12, ylabelsize=12,
        # xticks=LogTicks(WilkinsonTicks(8)),
        xticks = LogTicks(-9:3:0), yticks = LogTicks(-15:5:0),
        xscale=log10, yscale=log10, yminorticksvisible=false,
        xticklabelsize=9, yticklabelsize=9
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
    db = filter(row -> isfile(freqdatadirname*"frequencydata_$(row.environmentname).csv"), db)

    #/ Def. Taylor's law
    σf(μ, a) = @. a[1] * μ^a[2]
    #~ Allocate 
    afit = Array{Float64}(undef, length(db.environmentname))
    γfit = Array{Float64}(undef, length(db.environmentname))

    for (i, envname) in enumerate(db.environmentname)
        #/ Load relative frequency data
        fname = freqdatadirname * "/frequencydata_$(envname).csv"
        edb = CSV.read(fname, DataFrame, delim=", ")

        if !raw
            #/ Do some operations to 'bin' the frequencies as not to crowd the plot
            #  with individual scatter points
            #~ Bin using the `CategoricalArrays.cut`-function
            logmin, logmax = extrema(edb[!,:log_frequency])
            bins = range(logmin, logmax, nbins)
            edb[!,:bin] = cut(edb[!,:log_frequency], bins, extend=true)
            #~ Compute the mean of the variance for each of the bins
            bdb = @chain edb begin
                @by(:bin, :mean_variance = Statistics.mean(:var_frequency))
            end
            #~ Compute x and y for plotting
            μplot = exp.(collect(bins)[1:end-1] .+ diff(bins))
            σplot = bdb[!,:mean_variance]
            
            #/ Estimate params Taylor's law, and store for later
            afit[i], γfit[i] = fit_taylor(edb)
            
            #/ Plot
            scatter!(
                ax, μplot, σplot, markersize=4, strokewidth=.5,
                color=colors[i], marker=markers[i] #, label=envname
            )
        else
            #/ Simply plot the raw means and variances
            μplot = edb[!,:mean_frequency]
            σplot = edb[!,:var_frequency]
            #/ Plot
            scatter!(
                ax, μplot, σplot, markersize=3, strokewidth=.5,
                color=colors[i], marker=markers[i], label=envname
            )
        end
    end

    #/ Plot (fitted) Taylor's law
    xplot = exp10.(-9:3:0)
    amean = mean(afit)
    γmean = mean(γfit)
    powerlaw = lines!(
        ax, xplot, σf(xplot,[amean,γmean]), linewidth=.8, color=:black
    )
    label = L"\sigma^2 \propto \mu^{%$(round(γmean,digits=2))}"
    #~ Add legend
    axislegend(
        ax, [powerlaw], [label], position=:lt, framevisible=false, labelsize=11,
        patchsize=(5,1), padding=0, margin=(3,0,0,3), patchlabelgap=3,
    )

    #/ Resize and save (if desired)
    resize_to_layout!(fig)
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

########################
### HELPER FUNCTIONS ###
"Find the parameters of Taylor's law given a DataFrame"
function fit_taylor(db::DataFrame)
    x = db[!,:mean_frequency]
    y = db[!,:var_frequency]
    sol = CurveFit.power_fit(x, y)
    return sol
end

"Find the parameters of Taylor's law given x and y data"
function fit_taylor(x::Array{Float64}, y::Array{Float64})
    sol = CurveFit.power_fit(x, y)
    return sol
end

end # module MADPlotter
#/ End module
