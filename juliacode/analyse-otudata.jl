#= Simple module that analyses the raw OTU data
   Contains functions for loading and saving data, and relies on other modules to perform
   the more specialized computations, such as computing the moments.

#!note; assumes a directory tree as;
├── data
   ├── csv
   ├── jld
   └── rdata
├── juliacode
   ├── analyse-otudata.jl
   ├── compute-histogram.jl
   └── compute-moments.jl

where `data/` will contain both the raw and the filtered/analysed data. When the tree differs
significantly, the `const` paths should be changed accordingly.

#!note; for plotting the data, see `plot-mad.jl`
=#
#/ Start module
module OTUData

#/ Packages
using CSV, JLD2
using RData
using Chain, DataFrames, DataFramesMeta
using Distributions, Statistics

#/ Local packages
include("compute-histogram.jl")
include("compute-moments.jl")
using .Histogram
using .Moments

#~ Specify paths
#!note: if these do not exist, create them 
const RDATAPATH = "../data/rdata/"
const CSVDATAPATH = "../data/csv/crosssectional/"
const JLDATAPATH = "../data/jld/crosssectional/"
map(mkpath, [RDATAPATH, CSVDATAPATH, JLDATAPATH])

#################
### FUNCTIONS ###
"Load, split and filter data, and afterwards compute statistics for each environment"
function analyse(;
    rdatafilename = RDATAPATH*"crosssecdata.RData",
    split=true,     #~ Flag to split raw data into environment-specific data
    filter=true,    #~ Flag to filter raw data based on counts, reads, etc.
    compute=true,   #~ Flag to compute statistics (mean, var, etc.) from (filtered) data
    sad=true,       #~ Flag to compute histogram of abundances of (filtered) data
    mad=true,       #~ Flag to compute histogram of mean abundances of (filtered) data
    afd=true,       #~ Flag to compute histogram of abundance fluctuations of (filtered) data
    moments=true,   #~ Flag to compute estimates of moments of (filtered) data
    dry=false       #~ Flag for a 'dry' run, wherein nothing is saved (may break things)
)
    #/ Load and split data
    if split 
        envnamesdb = split_data(load_data(; rdatafilename=rdatafilename), dry=dry)
    else 
        envnamesdb = CSV.read(
            CSVDATAPATH*"environmentnames.csv", DataFrame, delim=", ",
            types=Dict(:projectclassification => String, :environmentname => String)
        )
    end
    #/ Load cutoffs
    cutoffsdb = CSV.read(
        CSVDATAPATH*"cutoffs.csv", DataFrame, delim=", ",
        types=Dict(:environmentname => String, :cutoff => Float64)
    )
    
    #/ Loop through all environments, filter the data, compute moments, and store
    for (i,env) in enumerate(envnamesdb.environmentname)
        if isfile(CSVDATAPATH * "rawotudata_$(env).csv")
            @info "Analysing... [env: $(env)]"
            edb = CSV.read(CSVDATAPATH * "rawotudata_$(env).csv", DataFrame, delim=", ")
            if filter
                #/ Filter 
                edb = filter_data(edb)
            end
            #!note: if edb=nothing, then there are no datapoints that 'survive' the filtering
            if !isnothing(edb)
                #/ if compute=true, compute summary statistics
                #!note: Using the (filtered) statsdb both the lognormal and Taylor's law
                #       can be extracted, so return the DataFrame here for completeness
                if compute
                    #/ Compute statistics
                    cutoff = only(@subset(cutoffsdb, :environmentname .== env)[!,:cutoff])
                    statsdb = compute_summarystatistics(edb, cutoff=cutoff)
                    logfreqdb = compute_logfrequencies(edb, cutoff=cutoff)
                    rescaledlogfreqdb = compute_rescaledlogfrequencies(edb, cutoff=cutoff)
                    # return statsdb, logfreqdb, rescaledlogfreqdb
                    #/ Save
                    if !dry
                        statsfname = CSVDATAPATH*"meanfrequencydata_$(env).csv"
                        CSV.write(statsfname, statsdb, delim=", ")
                        logfreqfname = CSVDATAPATH * "logfrequencydata_$(env).csv"
                        rescaledlogfreqfname = CSVDATAPATH *
                                               "rescaledlogfrequencydata_$(env).csv"
                        #~ write unscaled frequencies
                        __logfreqdb =
                            DataFrames.select(logfreqdb, [:otu_id,:log_frequency])
                        CSV.write(logfreqfname, __logfreqdb, delim=", ")
                        #~ write rescaled frequencies
                        __rescaledlogfreqdb =
                            DataFrames.select(rescaledlogfreqdb, [:otu_id,:log_frequency])
                        CSV.write(rescaledlogfreqfname, __rescaledlogfreqdb, delim=", ")
                    end
                else
                    statsfname = CSVDATAPATH*"meanfrequencydata_$(env).csv"
                    statsdb = CSV.read(statsfname, DataFrame, delim=", ")
                    logfreqfname = CSVDATAPATH*"logfrequencydata_$(env).csv"
                    rescaledlogfreqfname = CSVDATAPATH*"rescaledlogfrequencydata_$(env).csv"
                    rescaledlogfreqdb = CSV.read(rescaledlogfreqfname, DataFrame, delim=", ")
                end
                #/ if sad=true, compute the histogram of the abundance distribution
                #  without taking into account the mean across samples
                #!note: the relevant database is
                #!note: the relevant column is `:log_frequency`
                if sad
                    fh = Histogram.compute_fhist(logfreqdb[!,:log_frequency])
                    if !dry
                        JLD2.jldsave(JLDATAPATH * "sadfhist_$(env).jld2"; histogram=fh)
                    end
                end
                #@TODO IMPLEMENT THIS
                #/ if mad=true, compute the histogram of mean log frequencies
                #!note: the relevant column is `:mean_log_frequency`
                if mad
                    #/ Compute histogram
                    fh = Histogram.compute_fhist(statsdb[!,:mean_log_frequency])
                    if !dry
                        JLD2.jldsave(JLDATAPATH * "madfhist_$(env).jld2"; histogram = fh)
                    end
                end
                #/ if afd=true, compute the histogram of log frequencies accross samples
                #!Note: the relevant column is `:log_frequency`
                if afd
                    #/ Compute histogram
                    fh = Histogram.compute_fhist(rescaledlogfreqdb[!,:log_frequency])
                    if !dry
                        JLD2.jldsave(JLDATAPATH * "afdfhist_$(env).jld2"; histogram = fh)
                    end
                end
            end
        else
            @info "No data available. [env: $(env)]"
        end
    end

    #/ If the moments are to be computed, loop again and compute for all environments
    #  the moments using the desired procedure(s)
    if moments
        statsfname = CSVDATAPATH * "environmentstats.csv"
        if !isfile(statsfname)
            envstatsdb = @transform(envnamesdb, :mu=0.0, :sigma=0.0, :cutoff=0.0)
        else
            envstatsdb = CSV.read(statsfname, DataFrame, delim=", ")
        end
        for (i,env) in enumerate(envnamesdb.environmentname)
            cutoff = only(@subset(cutoffsdb, :environmentname .== env)[!,:cutoff])
            freqdatafname = CSVDATAPATH*"meanfrequencydata_$(env).csv"
            if isfile(freqdatafname)
                try 
                    _idx = findall(envstatsdb.environmentname.==env)[begin]
                    db = CSV.read(freqdatafname, DataFrame, delim=", ")
                    #/ Compute moments
                    logfreqs = db[!,:mean_log_frequency]
                    μest, σest = Moments.fittrunclognormal(
                        logfreqs;
                        uguess = [minimum(logfreqs), std(logfreqs), cutoff],
                        lower = cutoff
                    )
                    envstatsdb[_idx,:mu] = μest
                    envstatsdb[_idx,:sigma] = σest
                    envstatsdb[_idx,:cutoff] = cutoff
                catch e
                    @info "error" e
                    println("frequency data exists but datafile is empty, skipping [$(env)]")
                    continue
                end
            end
        end
        envstatsdb = @subset(envstatsdb, :sigma .> 0.0, :cutoff .> -Inf)
        CSV.write(statsfname, envstatsdb, delim=", ")
    end    
    @info "Done."
end

########################
### HELPER FUNCTIONS ###
"Load the RData as a DataFrame"
function load_data(; rdatafilename = RDATAPATH*"crosssecdata.RData")
    @info "Loading raw data..."
    db = RData.load(rdatafilename)["datatax"]
    db = @transform(db, :classification = String.(:classification))
    return db
end

"""
   Split data into data for distinct environments, such that for each environment we
   have a separate database that can be analyzed
"""
function split_data(db::DataFrame; dry=false)
    @info "Splitting raw data..."
    #/ 1. Create unique ID for project and classification
    db = @transform(db, :projectclassification = :project_id .* :classification)
    #~ Load short-hand names for the environments
    #  these are defined in CSVDATAPATH/environmentnames.csv
    environmentnamesdb = CSV.read(
        CSVDATAPATH*"environmentnames.csv", DataFrame, delim=", ",
        types=Dict(:projectclassification => String, :environmentname => String)
    )
    #~ Replace the :project_id and :classification columns by a single column
    #  this makes it easier down the line, and also provides the :environmentname column
    db = @chain begin
        innerjoin(db, environmentnamesdb, on=:projectclassification)
        select(Not([:project_id, :classification, :projectclassification]))
    end

    #/ 2. For each environmentname, select the subset of data from that environment and
    #     store it seperately for further analysis
    for environment in environmentnamesdb.environmentname
        edb = @subset(db, :environmentname .== environment)
        #~ Save dataframe
        if !dry 
            filename = CSVDATAPATH * "rawotudata_$(environment).csv"
            CSV.write(filename, edb, delim=", ")
        end
    end
    return environmentnamesdb
end

function filter_data(db::DataFrame;
    minsamples = 30,
    minreads = 10_000,
    mincounts = 1,
    remove_runs = ["ERR1104477", "ERR1101508", "SRR2240575"] # bad runs filtered by Grilli
)
    #/ Check the total number of samples
    #~ if not enough samples, do nothing (i.e., skip it)
    nsamples = length(unique(db[!,:sample_id]))
    if nsamples < minsamples
        @info "not enough samples, skipping [env: $(first(db[!,:environmentname]))]"
        return nothing
    end

    #/ Filter entries and create filtered dataframe
    fdb = @chain db begin
        @rsubset(!in(:run_id, remove_runs))
        @subset(:nreads .> minreads)
        @subset(:count .> mincounts)
    end
    #/ Return filtered dataframe only when non-empty
    (nrow(fdb) > 0) && (return fdb)
    #~ Otherwise, return nothing
    @info "not enough reads or counts [env: $(first(db[!,:environmentname]))]"
    return nothing
end

"""
Compute log frequencies (relative abundances) across communities (samples)

!note: Simply aggregates all frequencies into a single distribution without rescaling       
!note: Assumes that the frequency (relative abundance) is defined as #count/#totalreads
"""
function compute_logfrequencies(fdb::DataFrame; cutoff = -100.0)
	  #/ Compute the total number of runs/samples for the specific environment
    nruns = length(unique(fdb[!,:run_id]))
    #/ Chain multiple dataframe operations to compute the rescaled log frequency
    #~ rescaled log frequency = log((x - μ)/σ)
    db = @chain fdb begin
        #~ Compute frequencies
        @transform(:frequency = :count ./ :nreads)
        @transform(:log_frequency = log.(:frequency))
        @subset(:log_frequency .> cutoff)
        @select(:otu_id, :sample_id, :run_id, :log_frequency)
    end

    db = @chain db begin
        #~ Omit (log) frequencies that are NaN and/or missing
        @transform(:log_frequency = coalesce.(:log_frequency, -Inf))
        @subset(:log_frequency .> -Inf, :log_frequency .< Inf)
        @subset(:log_frequency .> cutoff)
    end

    return db
end

"""
Compute statistics for a specific filtered dataframe

!note: Assumes that the frequency is defined as #count/#totalreads
"""
function compute_summarystatistics(fdb::DataFrame; cutoff::Float64 = -100.0)
    #/ Compute the total number of runs/samples for the specific environment
    nruns = length(unique(fdb[!,:run_id]))
    #/ Chain multiple dataframe operations to compute mean frequency
    statsdb = @chain fdb begin
        #~ Compute frequencies
        @transform(:frequency = :count ./ :nreads)
        #~ For each species (OTU), group, and compute mean and variance of the frequency
        @by(
            :otu_id,
            :mean_frequency = Statistics.mean(skipmissing(:frequency)),
            :var_frequency = Statistics.var(skipmissing(:frequency), corrected=false),
            #~ Compute the occupancy; i.e. the fraction of samples (runs) wherein the focal
            #  species is actually present
            :occupancy = length(:otu_id) ./ nruns
        )
        #~ Take the occupation number into account
        #~ this means that μ → o⋅μ and σ² → o⋅[σ²+μ²(1-o)], where o the occupancy
        @transform(:mean_frequency = :mean_frequency .* :occupancy)
        @transform(:var_frequency = :var_frequency .+ :mean_frequency.^2 .* (1 .- :occupancy))
        @transform(:var_frequency = :var_frequency .* :occupancy)
        #~ Perform a log-transform on the mean-frequency (needed for lognormal)
        @transform(:mean_log_frequency = log.(:mean_frequency))
        #~ Select only those above a specified cutoff (filter)
        @subset(:mean_log_frequency .> cutoff, :var_frequency .> 0.0)
    end
    return statsdb
end

"""
Compute rescaled log frequencies (relative abundances) across communities (samples)

!note: Assumes that the frequency is defined as #count/#totalreads
"""
function compute_rescaledlogfrequencies(fdb::DataFrame; cutoff = -100.0)
    #/ Compute the total number of runs/samples for the specific environment
    nruns = length(unique(fdb[!,:run_id]))
    #/ Chain multiple dataframe operations to compute the rescaled log frequency
    #~ rescaled log frequency = log((x - μ)/σ)
    db = @chain fdb begin
        #~ Compute frequencies
        @transform(:frequency = :count ./ :nreads)
        @transform(:log_frequency = log.(:frequency))
        @subset(:log_frequency .> cutoff)
    end

    #/ Compute summary statistics for each otu_id
    summarydb = @chain db begin
        @by(
            :otu_id,
            :mean_logfrequency = mean(skipmissing(:log_frequency)),
            :std_logfrequency = std(skipmissing(:log_frequency), corrected=false),
            :occupancy = length(:otu_id) ./ nruns
        )
        @subset(:std_logfrequency .> 0.0, :occupancy .≈ 1.0)
        @select(:otu_id, :mean_logfrequency, :std_logfrequency)
    end

    #/ Rescale the log frequency by the summary statistics
    #!note: `missing` values are propagated and need to be filtered out
    db = DataFrames.leftjoin(db, summarydb, on=:otu_id)
    db = @chain db begin
        @transform(:log_frequency = (:log_frequency.-:mean_logfrequency)./:std_logfrequency) 
        #~ Omit (log) frequencies that are NaN and/or missing
        @transform(:log_frequency = coalesce.(:log_frequency, -Inf))
        @subset(:log_frequency .> -Inf, :log_frequency .< Inf)
        @subset(:log_frequency .> cutoff)
    end
    return db
end

end # module OTUData
#/ End module
