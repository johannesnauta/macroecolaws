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
module lOTUData

#/ Packages
using CSV, JLD2
using RData
using Chain, DataFrames, DataFramesMeta
using Distributions, Statistics

#/ Local packages
include("compute-histogram.jl")
include("compute-moments.jl")
include("compute-mixture.jl")
using .Histogram
using .Moments
using .Mixture

#~ Specify paths
#!note: if these do not exist, create them 
const RDATAPATH = "../data/rdata/"
const CSVDATAPATH = "../data/csv/longitudinal/"
const JLDATAPATH = "../data/jld/longitudinal/"
map(mkpath, [RDATAPATH, CSVDATAPATH, JLDATAPATH])

#################
### FUNCTIONS ###
"Load, split and filter data, and afterwards compute statistics for each environment"
function analyse(;
    rdatafilename = RDATAPATH*"longitudinal.RData",
    split=true,      #~ Flag to split raw data into environment-specific data
    filter=true,     #~ Flag to filter raw data based on counts, reads, etc.
    compute=true,    #~ Flag to compute statistics (mean, var, etc.) from (filtered) data
    mad=true,        #~ Flag to compute histogram of mean abundances of (filtered) data
    afd=true,        #~ Flag to compute histogram of abundance fluctuations of (filtered) data
    moments=true,    #~ Flag to compute estimates of moments of (filtered) data
    fitparams=true,  #~ Flag to individual parameters for the (log) frequency distribution(s)
    pearson=true,    #~ Flag to compute Pearsons correlation coefficient
    dry=false        #~ Flag for a 'dry' run, wherein nothing is saved [note: may break things]
)
    #/ Load and split data
    if split
        # return load_data(; rdatafilename=rdatafilename)
        envnamesdb = split_data(load_data(; rdatafilename=rdatafilename), dry=dry)
    else 
        envnamesdb = CSV.read(
            CSVDATAPATH*"environmentnames.csv", DataFrame, delim=", ",
            types=Dict(:host_id => String, :samplesite => String, :environmentname => String)
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
                    # return compute_logfrequencies(edb, cutoff=cutoff)
                    statsdb = compute_summarystatistics(edb, cutoff=cutoff)
                    logfreqdb = compute_logfrequencies(edb, cutoff=cutoff)
                    rescaledlogfreqdb = compute_rescaledlogfrequencies(edb, cutoff=cutoff)
                    #/ Save
                    if !dry
                        statsfname = CSVDATAPATH*"meanfrequencydata_$(env).csv"
                        CSV.write(statsfname, statsdb, delim=", ")
                        logfreqfname = CSVDATAPATH * "logfrequencydata_$(env).csv"
                        rescaledlogfreqfname = CSVDATAPATH *
                                               "rescaledlogfrequencydata_$(env).csv"
                        #~ write unscaled frequencies
                        __logfreqdb = DataFrames.select(
                            logfreqdb, [:otu_id,:experiment_day,:log_frequency]
                        )
                        CSV.write(logfreqfname, __logfreqdb, delim=", ")
                        #~ write rescaled frequencies
                        __rescaledlogfreqdb = DataFrames.select(rescaledlogfreqdb,
                            [:otu_id,:experiment_day,:log_frequency,:sample_log_frequency]
                        )
                        CSV.write(rescaledlogfreqfname, __rescaledlogfreqdb, delim=", ")
                    end
                else
                    statsfname = CSVDATAPATH*"meanfrequencydata_$(env).csv"
                    statsdb = CSV.read(statsfname, DataFrame, delim=", ")
                    logfreqfname = CSVDATAPATH*"logfrequencydata_$(env).csv"
                    rescaledlogfreqfname = CSVDATAPATH*"rescaledlogfrequencydata_$(env).csv"
                    rescaledlogfreqdb = CSV.read(rescaledlogfreqfname, DataFrame, delim=", ")
                end
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
                    sfh = Histogram.compute_fhist(rescaledlogfreqdb[!,:sample_log_frequency])
                    if !dry
                        JLD2.jldsave(JLDATAPATH * "afdfhist_$(env).jld2"; histogram=fh)
                        JLD2.jldsave(JLDATAPATH * "sampleafdfhist_$(env).jld2"; histogram=sfh)
                    end
                end
                #/ if pearson=true, compute the Pearson correlation coefficient between each
                #  of the OTUs that have enough (overlapping) datapoints
                if pearson
                    pearsondb = compute_pearson(logfreqdb)
                    if !dry
                        JLD2.jldsave(
                            JLDATAPATH * "pearsoncorrelation_$(env).jld2"; pearson=pearsondb
                        )
                    end
                end
                #/ if fitparams=true, compute the parameters of the distribution for each
                #  of the OTUs that have enough datapoints
                #!note: assumes that the distribution is a Gamma distribution by default
                if fitparams
                    fitdb = compute_shapescale(logfreqdb)
                    distfitfname = CSVDATAPATH * "distfitdata_$(env).csv"
                    CSV.write(distfitfname, fitdb, delim=", ")
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
                    # return db
                    #/ Compute moments
                    logfreqs = db[!,:mean_log_frequency]
                    μest, logσest = Moments.fittrunclognormal(
                        logfreqs; uguess = [mean(logfreqs),log(std(logfreqs))], lower=cutoff
                    )
                    envstatsdb[_idx,:mu] = μest
                    envstatsdb[_idx,:sigma] = exp(logσest)
                    envstatsdb[_idx,:cutoff] = cutoff
                catch e
                    println(e)
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
    db = RData.load(rdatafilename)["proj_time"]
    db = @transform(
        db,
        :classification = String.(:classification),
        :host_id = String.(:host_id),
        :samplesite = String.(:samplesite),
        :date_collection = String.(:date_collection)
    )
    return db
end

"""
   Split data into data for distinct environments, such that for each environment we
   have a separate database that can be analyzed
"""
function split_data(db::DataFrame; dry=false)
    @info "Splitting raw data..."
    #/ 1. Create unique ID for host ID and sampling site
    db = @chain db begin
        @transform(:site_id = :host_id .* :samplesite)
        @select(Not([:project_id, :classification]))
    end
    #~ Load short-hand names for the environments
    #  these are defined in CSVDATAPATH/environmentnames.csv
    environmentnamesdb = CSV.read(
        CSVDATAPATH*"environmentnames.csv", DataFrame, delim=", ",
        types=Dict(:host_id => String, :samplesite => String, :environmentname => String)
    )
    environmentnamesdb = @chain environmentnamesdb begin
        @transform(:site_id = :host_id .* :samplesite)
        @select(Not([:host_id, :samplesite]))
    end
    
    #~ Replace the :project_id and :classification columns by a single column
    #  this makes it easier down the line, and also provides the :environmentname column
    db = innerjoin(db, environmentnamesdb, on=:site_id)
    db = select(db, Not([:host_id, :samplesite]))

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

"""
    Filter the data by omitting all samples/environments that fall below some thresholds
"""
function filter_data(db::DataFrame;
    minsamples = 30,
    minreads = 1000,
    mincounts = 1
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
function compute_summarystatistics(
    fdb::DataFrame;
    cutoff::Float64 = -100.0,
    minoccurrences::Int = 1
)
    #/ Compute the total number of runs/samples for the specific environment
    nruns = length(unique(fdb[!,:run_id]))
    #/ Chain multiple dataframe operations to compute mean frequency
    statsdb = @chain fdb begin
        #~ Compute frequencies
        @transform(:frequency = :count ./ :nreads)
        @transform(:log_frequency = log.(:frequency))
        #~ Filter frequencies that are too rate
        @subset(:log_frequency .> cutoff)
        #~ For each species (OTU), group, and compute mean and variance of the frequency
        @by(
            :otu_id,
            :noccurrences = length(:otu_id),
            :mean_frequency = Statistics.mean(skipmissing(:frequency)),
            :var_frequency = Statistics.var(skipmissing(:frequency), corrected=false),
            #~ Compute the occupancy; i.e. the fraction of samples (runs) wherein the focal
            #  species is actually present
            :occupancy = length(:otu_id) ./ nruns
        )
        #~ Take the occupation number into account
        #~ i.e., only take into account species that appeared in most runs
        @subset(:noccurrences .> minoccurrences, :var_frequency .> 0.0)
        #~ this means that μ → o⋅μ and σ² → o⋅[σ²+μ²(1-o)], where o the occupancy
        @transform(:mean_frequency = :mean_frequency .* :occupancy)
        @transform(:var_frequency = :var_frequency .+
                                        :mean_frequency.^2 .* (1 .- :occupancy))
        @transform(:var_frequency = :var_frequency .* :occupancy)
        @transform(:mean_log_frequency = log.(:mean_frequency))
        @transform(:std_frequency = sqrt.(:var_frequency))
        @select(Not(:var_frequency))
        #~ Perform a log-transform on the mean-frequency (needed for lognormal)
        # @transform(:mean_log_frequency = log.(:mean_frequency))
        #~ Select only those with non-zero variance (i.e., >1 reads)
        # @subset(:std_frequency .> 0.0)
    end
    return statsdb
end

"""
Compute log frequencies (relative abundances) across communities (samples)

!note: Simply aggregates all frequencies into a single distribution without rescaling       
!note: Assumes that the frequency (relative abundance) is defined as #count/#totalreads
"""
function compute_logfrequencies(fdb::DataFrame; cutoff = -100.0)
    #/ Chain multiple dataframe operations to compute the log frequency
    db = @chain fdb begin
        #~ Compute frequencies
        @transform(:frequency = :count ./ :nreads)
        @transform(:log_frequency = log.(:frequency))
        #~ ensure all log frequencies are above a cutoff
        @subset(:log_frequency .> cutoff)
        @select(:otu_id, :sample_id, :run_id, :experiment_day, :log_frequency)
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
Compute rescaled log frequencies (relative abundances) across communities (samples)

!note: Assumes that the frequency is defined as #count/#totalreads
"""
function compute_rescaledlogfrequencies(
    fdb::DataFrame;
    cutoff = -100.0,
    minoccurrences::Int = 50
)
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

    #~ filter out entries that do not have enough occurences to compute a distribution
    #~ and those that have many repeated values (i.e., where the min=max)
    odb = @by(
        db, :otu_id,
        :noccurrences = length(:otu_id),
        :min_log_frequency = minimum(:log_frequency),
        :max_log_frequency = maximum(:log_frequency)
    )
    db = leftjoin(db, odb, on=:otu_id)
    db = @subset(db,
        :noccurrences .> minoccurrences,
        :max_log_frequency .- :min_log_frequency .> 0
    )

    #/ Compute summary statistics for each otu_id
    summarydb = @chain db begin
        @by(
            :otu_id,
            :mean_log_frequency = mean(Histogram.compute_fhist(:log_frequency)),
            :var_log_frequency = std(Histogram.compute_fhist(:log_frequency)).^2,
            :sample_mean_log_frequency = mean(:log_frequency),
            :sample_var_log_frequency = std(:log_frequency, corrected=false).^2,
            :max_log_frequency = maximum(:log_frequency),
            :min_log_frequency = minimum(:log_frequency),
            :occupancy = length(:otu_id) ./ nruns
        )
        @subset(:var_log_frequency .> 0.0, :sample_var_log_frequency .> 0.0)
        #~ Take the occupancy into account
        @subset(:occupancy .≈ 1.0)
        # @transform(:mean_log_frequency = :mean_log_frequency .* :occupancy)
        # @transform(:var_log_frequency = :var_log_frequency .+
        #                                 :mean_log_frequency.^2 .* (1 .- :occupancy))
        # @transform(:var_log_frequency = :var_log_frequency .* :occupancy)
        # @transform(:sample_mean_log_frequency = :sample_mean_log_frequency .* :occupancy)
        # @transform(:sample_var_log_frequency = :sample_var_log_frequency .+
        #                                 :sample_mean_log_frequency.^2 .* (1 .- :occupancy))
        # @transform(:sample_var_log_frequency = :sample_var_log_frequency .* :occupancy)
        #~ Compute the standard deviation
        @transform(:std_log_frequency = sqrt.(:var_log_frequency))
        @transform(:sample_std_log_frequency = sqrt.(:sample_var_log_frequency))
        #~ Select only relevant columns
        @select(
            :otu_id,
            :mean_log_frequency, :std_log_frequency,
            :sample_mean_log_frequency, :sample_std_log_frequency
        )
    end
    
    #/ Rescale the log frequency by the summary statistics
    #~ @TODO Is this needed for the AFD??
    #!note: `missing` values are propagated and need to be filtered out
    db = DataFrames.leftjoin(db, summarydb, on=:otu_id)
    db = @chain db begin
        #~?Why does this particular step 'work', as the log_frequencies are *not*
        #  normally distributed (in fact, they are most likely gamma distributed), so the
        #  mean here has no statistical meaning
        @transform(
            :sample_log_frequency =
                (:log_frequency.-:sample_mean_log_frequency) ./ :sample_std_log_frequency
        )
        @transform(:log_frequency = (:log_frequency.-:mean_log_frequency)./:std_log_frequency)
        #~ Omit (log) frequencies that are NaN and/or missing
        @transform(:log_frequency = coalesce.(:log_frequency, -Inf))
        @subset(:log_frequency .> -Inf, :log_frequency .< Inf)
        @transform(:sample_log_frequency = coalesce.(:sample_log_frequency, -Inf))
        @subset(:sample_log_frequency .> -Inf, :sample_log_frequency .< Inf)
        @subset(:log_frequency .> cutoff, :sample_log_frequency .> cutoff)
    end
    return db
end

"""
Compute shape and scale parameters of OTUs with at least a specific number of days sampled
"""
function compute_shapescale(fdb::DataFrame; mindays::Int = 30, prior=Distributions.Gamma)
    #~ Compute the total no. of days that the OTU was measured
    daydb = @chain fdb begin
        @by(:otu_id, :ndays = length(unique(:experiment_day)))
        @subset(:ndays .> mindays)
    end
    fdb = rightjoin(fdb, daydb, on=:otu_id)

    #~ For each otu_id, fit a Gamma mixture distribution and
    #  collect the shape and scale parameters
    pdb = @chain fdb begin
        @transform(:frequency = exp.(:log_frequency))
        @groupby(:otu_id)
        @combine(
            :mleparams = Mixture.fit_mixture(
                :frequency, prior=prior,
                guess = [mean(:frequency) / 1., 1., 0.9, maximum(:frequency), 1.]
            )
        )
        #~ Get the field of the :mleparams struct with the symbol ^(:symbol)
        #  otherwise DataFramesMeta macros thing :symbol is a column, which they are not
        @transform(:shape = getfield.(:mleparams, ^(:shape)))
        @transform(:scale = getfield.(:mleparams, ^(:scale)))
        @transform(:ε = getfield.(:mleparams, ^(:ε)))
        @select(Not(:mleparams))
    end
    return pdb
end

"""
Compute Pearson correlation coefficient between trajectories
"""
function compute_pearson(fdb::DataFrame; mindays::Int = 50)
    #~ Compute the total no. of days that the OTU was measured
    daydb = @chain fdb begin        
        @by(:otu_id, :ndays = length(unique(:experiment_day)))
        @subset(:ndays .> mindays)
    end
    fdb = rightjoin(fdb, daydb, on=:otu_id)

    #/ For each otu_id, compute the Pearson correlation coefficient with all other OTUs
    #~ compute statistics, mean and variance, in order to compute the correlation
    summarydb = @chain fdb begin
        @by(
            :otu_id,
            :sample_mean_log_frequency = mean(:log_frequency),
            :sample_std_log_frequency = std(:log_frequency, corrected=false)
        )
    end

    #/ Extract all timeseries
    seriesdb = @chain fdb begin
        @select(:experiment_day, :otu_id, :log_frequency)
        unstack(:experiment_day, :otu_id, :log_frequency)
    end

    #/ Allocate
    pearsondb = DataFrame(otu_id1 = Int[], otu_id2 = Int[], pearson=Float64[])

    #/ Compute Pearson correlations
    otus = names(seriesdb)[2:end]
    for i in eachindex(otus), j in (i+1):length(otus)        
        #~ note: there is most likely a DataFrame-y way to do things, but this if more clear
        x = seriesdb[!,otus[i]]
        y = seriesdb[!,otus[j]]
        iidx = findfirst(==(parse(Int,otus[i])), summarydb.otu_id)
        jidx = findfirst(==(parse(Int,otus[j])), summarydb.otu_id)        
        meanx = summarydb[iidx,:].sample_mean_log_frequency
        stdx = summarydb[iidx,:].sample_std_log_frequency
        meany = summarydb[jidx,:].sample_mean_log_frequency
        stdy = summarydb[jidx,:].sample_std_log_frequency
        #~ compute Pearson correlation
        xshift = passmissing(/).(passmissing(-).(x, meanx), stdx)
        yshift = passmissing(/).(passmissing(-).(y, meany), stdy)
        s = passmissing(*).(xshift, yshift)
        n = count(!ismissing, s)
        if n > mindays
            push!(
                pearsondb,
                (
                    otu_id1 = parse(Int,otus[i]),
                    otu_id2 = parse(Int,otus[j]),
                    pearson = sum(skipmissing(s)) / (n - 1)
                )
            )
        end
    end
    #~ Filter out entries that have infinite correlation as n=1, or when |ρ|>1
    #~ note: should not be needed if mindays is large enough, but do it anyways
    pearsondb = @subset(pearsondb, :pearson .< 1., :pearson .> -1.)
    
    return pearsondb
end

end # module OTUData
#/ End module
