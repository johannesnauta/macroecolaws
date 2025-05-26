#= Module for plotting mean abundance distributions (MADs) from data

!note: some functionality is limited (such as no diff. markers), as the code was originally 
       written using a specific personalized theme that is not included in this repository
=#
#/ Start module
module lPlotter

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
		Plot abundance trajectories of the first n OTUs that have the most datapoints
"""
function plot_trajectories(;
    topno::Int = 1,
    prefix::String = "longitudinal/",
    trajectorydir::String = CSVDATAPATH * prefix,
    envnamefname::String = CSVDATAPATH * prefix * "environmentnames.csv"
)
    #/ Create figure
    width = 1.5 * 246
    fig = Figure(;
        size=(width,width/1.87), figure_padding=(2,8,2,2), backgroundcolor=:transparent
    )

    #/ Load environmentnames
    envdb = CSV.read(envnamefname, DataFrame, delim=", ", types=String)
    #/ Filter environments to include only those for which trajectories exist
    envdb = filter(
        row -> isfile(trajectorydir*"logfrequencydata_$(row.environmentname).csv"), envdb
    )
    envdb = @transform(envdb, :texname = replace.(:environmentname, "_" => "\\;"))

    ax = [Axis(
        fig[i,j],
        title=L"\textrm{(%$(envdb.texname[j+4*(i-1)]))}", titlesize=7.5, titlegap=1,
        xlabel=L"\textrm{time (days)}", xlabelvisible=(i==2),
        ylabel=L"\textrm{log\;abundances}", ylabelvisible=(j==1),
        xlabelsize=8, ylabelsize=9,
        xticklabelsvisible=true, yticklabelsvisible=(j==1),
        xticks=LinearTicks(5),
        limits=(0,(i==1 ? 450 : 185),-11,0),
        xticklabelsize=7, yticklabelsize=7,
    ) for i in 1:2, j in 1:4]

    for (i, envname) in enumerate(envdb.environmentname)
        #/ Load trajectories
        db = CSV.read(trajectorydir * "logfrequencydata_$(envname).csv", DataFrame, delim=", ")
        #/ Select the top n trajectories with the most days
        #~ find the OTUs that have the most days
        topdb = @chain db begin
            @by(:otu_id, :ndays = length(unique(:experiment_day)))
            @orderby(-:ndays)
            first(topno)
        end
        (topno < Inf) && (topdb = first(topdb, topno))
        @info extrema(topdb.ndays)
        fdb = filter(row -> row.otu_id in topdb.otu_id, db)
        fdb = @orderby(fdb, :otu_id)

        #/ Plot
        for (n,otu) in enumerate(topdb.otu_id)
            x, y = div(i-1,4)+1, mod1(i,4)
            plotdb = filter(row -> row.otu_id == otu, fdb)
            plotdb = @orderby(plotdb, :experiment_day)
            lines!(
                ax[x,y], plotdb.experiment_day, plotdb.log_frequency, linewidth=.4
            )
        end
    end    

    return fig
end

"""
		Plot abundance trajectories of the first n OTUs that have the most datapoints, but this
    time only for a single environment
"""
function plot_trajectory(;
    topno::Int = 1,
    prefix::String = "longitudinal/",
    environment = "M_leftpalm",
    trajectorydir::String = CSVDATAPATH * prefix,
    envnamefname::String = CSVDATAPATH * prefix * "environmentnames.csv",
    augment = true, #~ augment with zeros when no data on a specific day
    savefig = true,
    figname = nothing
)
    #/ Create figure
    width = .9 * 246
    fig = Figure(;
        size=(width,width/1.67), figure_padding=(2,8,2,2), backgroundcolor=:transparent
    )

    #/ Load environmentnames
    envdb = CSV.read(envnamefname, DataFrame, delim=", ", types=String)
    #/ Filter environments to include only those for which trajectories exist
    envdb = filter(row -> row.environmentname == environment, envdb)
    envdb = @transform(envdb, :texname = replace.(:environmentname, "_" => "\\;"))

    ax = Axis(
        fig[1,1],
        title=L"\textrm{M\;left palm}", titlesize=10, titlegap=1,
        xlabel=L"\textrm{time (days)}",
        ylabel=L"\textrm{log\;rel.\;abundances}",
        xlabelsize=11, ylabelsize=11,
        xticklabelsvisible=true, yticklabelsvisible=true,
        xticks=LinearTicks(5),
        xminorticks=IntervalsBetween(5), xminorticksvisible=true, xminorticksize=1.5,
        yminorticks=IntervalsBetween(3), yminorticksvisible=true, yminorticksize=1.5,
        limits=(0,450,-9,0),
        xticklabelsize=8, yticklabelsize=8,
    )
    

    #/ Load trajectories
    db = CSV.read(trajectorydir * "logfrequencydata_$(environment).csv", DataFrame, delim=", ")
    #/ Select the top n trajectories with the most days
    #~ find the OTUs that have the most days
    topdb = @chain db begin
        @by(:otu_id, :ndays = length(unique(:experiment_day)))
        @orderby(-:ndays)
        first(topno)
    end
    (topno < Inf) && (topdb = first(topdb, topno))
    @info extrema(topdb.ndays)
    fdb = filter(row -> row.otu_id in topdb.otu_id, db)
    fdb = @orderby(fdb, :otu_id)

    #/ Plot
    for (n,otu) in enumerate(topdb.otu_id)
        plotdb = filter(row -> row.otu_id == otu, fdb)
        plotdb = @orderby(plotdb, :experiment_day)
        #~ augment with some zeros
        if augment
            start,stop = minimum(plotdb.experiment_day), maximum(plotdb.experiment_day)
            plotdays = vcat(plotdb.experiment_day, setdiff(start:stop,plotdb.experiment_day))
            sort!(plotdays)
            valdict = Dict(zip(plotdb.experiment_day, plotdb.log_frequency))
            plotfreqs = [get(valdict, day, -10) for day in plotdays]
            lines!(ax, plotdays, plotfreqs, linewidth=.3)
        else
            lines!(ax, plotdb.experiment_day, plotdb.log_frequency, linewidth=.3)
        end
    end
    (savefig && !isnothing(figname)) && (CairoMakie.save(figname, fig, pt_per_unit=1))
    return fig
end

end # module lPlotter
#/ End module
