# Macroecological laws in microbial systems
Repository to subtract macro-ecological laws from microbial abundance data.
At its core it aims to reproduce the results presented in [Grilli, 2020](https://www.nature.com/articles/s41467-020-18529-y). That is, it analyses the data to obtain
* a lognormal mean-abundance distribution (MAD)
* [not-yet-implemented] a gamma abundance fluctuation distribution (AFD)
* [not-yet-implemented] Taylor's law 

The codebase contains both the raw data and a submodule that is the original code that accompanied the paper.

## Running the `Julia` code

### Packages 
As this little repository is _not_ meant to be a functioning package/module, it does not include a `Project.toml` by which users can easily install the packages. 
The actual reason for this is that this code has been written as a part of a much larger codebase, and nesting `Project.toml`s within other (large) modules is a recipe for disaster. 
Therefore, it requires the user to install the appropriate packages themselves. 
There should be no dependence on specific versions (as far as I am aware).

The packages, in no particular order, that are used are
- I/O
    - `CSV.jl`, `JLD2.jl`, `RData.jl`
- dataframes
    - `Chain.jl`, `DataFrames.jl`, `DataFramesMeta.jl`
- statistics
    - `Statistics.jl`, `FHist.jl`
- solvers
    - `NonlinearSolve.jl`, `NaNMath.jl`, `SpecialFunctions.jl`, `Optim.jl`
- distributions and rng manipulation
    - `Distributions.jl`, `Random.jl`
- plotting
    - `CairoMakie.jl`, `LaTeXString.jl`

###
The main file to run is `juliacode/analyse-otudata.jl`, which (depending on the flags/options) does the following;
- [`split=true`] splits the data into separated OTUs for each of the "environments" (here, an "environment" is a specific projectclassification; names for the environments are defined in `data/csv/environmentnames.csv`)
- [`filter=true`] filters the split data by omitting entries which, for example, do not contain enough samples, counts, and or reads
- [`compute=true`] computes the the statistics of the data (e.g., log mean abundance, variance, etc.)
- [`histogram=true`] computes the histogram of the filtered data
- [`moments=true`] computes estimates of the moments (of the lognormal distribution), using either maximum-likelihood estimation (as done in the paper), or by fitting a (truncated) lognormal distribution
- [`dry=false`] if `true`, store the minimum no. of files (e.g., after splitting, nothing is stored, useful for debugging)

All paths are relative paths and the script assumes a directory tree as in this repository, which is
```bash 
├── data
	 ├── csv
   ├── jld
   └── rdata
├── juliacode
   ├── analyse-otudata.jl
   ├── compute-histogram.jl
   └── compute-moments.jl
```

### The `OTUData` module
The main script `analyse-otudata.jl` can be run by running in the `Julia` REPL (assumes `Revise.jl` is used)

```julia
julia> includet("analyse-otudata.jl")
julia> using .OTUData
julia> OTUData.analyse()
```

### Other modules
For plotting, see
- for plotting the lognormal distribution(s), see `plot/plot-mad.jl`
    - note that this module possibly relies on additional packages, such as `ElectronDisplay.jl` in order for the plots to be shown
- [...]


---
This small repository is meant to be shared in order to facilitate the recreation of the analysis of the OTU data provided by [Grilli, 2020](https://www.nature.com/articles/s41467-020-18529-y). 
The aim is to investigate how to recreate the macroecological laws, and to increase the likelihood of recreation by others who are, perhaps, less familiar with `R` and/or simply want something "that just runs". 
The intention of this script is therefore to be a standalone analyser and plotter for the specific OTU data provided and is by no means a full analysis. 
There may also be countless errors as, as of writing this, I have yet to reproduce the figures.
