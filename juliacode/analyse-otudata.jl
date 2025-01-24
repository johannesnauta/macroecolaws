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
using Statistics

#/ Local packages
include("compute-histogram.jl")
include("compute-moments.jl")
using .Histogram
using .Moments

#~ Specify paths
#!note: if these do not exist, create them 
const RDATAPATH = "../data/rdata/"
const CSVDATAPATH = "../data/csv/"
const JLDATAPATH = "../data/jld/"
map(mkpath, [RDATAPATH, CSVDATAPATH, JLDATAPATH])

#################
### FUNCTIONS ###
"Load, split and filter data, and afterwards compute statistics for each environment"
function analyse(;
    rdatafilename = RDATAPATH*"crosssecdata.RData",
    split=true,     #~ Flag to split raw data into environment-specific data
    filter=true,    #~ Flag to filter raw data based on counts, reads, etc.
    compute=true,   #~ Flag to compute statistics (mean, var, etc.) from (filtered) data
    histogram=true, #~ Flag to compute histogram of (filtered) data
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
                if compute
                    #/ Compute statistics
                    cutoff = only(@subset(cutoffsdb, :environmentname .== env)[!,:cutoff])
                    statsdb = compute_stats(edb, cutoff=cutoff)
                    #/ Save
                    if !dry
                        statsfname = CSVDATAPATH*"frequencydata_$(env).csv"
                        CSV.write(statsfname, statsdb, delim=", ")
                    end
                else
                    statsfname = CSVDATAPATH*"frequencydata_$(env).csv"
                    statsdb = CSV.read(statsfname, DataFrame, delim=", ")
                end
                if histogram
                    #/ Compute histogram
                    fh = Histogram.compute_fhist(statsdb[!,:log_frequency])
                    if !dry
                        JLD2.jldsave(JLDATAPATH * "fhist_$(env).jld2"; histogram = fh)
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
            freqdatafname = CSVDATAPATH*"frequencydata_$(env).csv"            
            if isfile(freqdatafname)
                db = CSV.read(freqdatafname, DataFrame, delim=", ")
                #/ Compute moments
                logfreqs = db[!,:log_frequency]
                μest, σest = Moments.fittrunclognormal(
                    logfreqs; uguess = [mean(logfreqs), std(logfreqs)], lower=cutoff
                )
                _idx = findall(envstatsdb.environmentname.==env)[begin]
                envstatsdb[_idx,:mu] = μest
                envstatsdb[_idx,:sigma] = σest
                envstatsdb[_idx,:cutoff] = cutoff
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
Compute statistics for a specific filtered dataframe

!note: Assumes that the frequency is defined as #count/#totalreads
"""
function compute_stats(fdb::DataFrame; cutoff::Float64 = -100.0)
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
            :num_samples = length(:sample_id),
            :occupation = length(:otu_id)      #~ similar to the n() function in `R`
        )
        #~ Take the occupation number into account
        @transform(:mean_frequency = :mean_frequency .* (:occupation ./ nruns))
        #~ Perform a log-transform on the mean-frequency (needed for lognormal)
        @transform(:log_frequency = log.(:mean_frequency))
        #~ Select only those above a specified cutoff (filter)
        @subset(:log_frequency .> cutoff)
    end
    return statsdb
end

end # module OTUData
#/ End module
