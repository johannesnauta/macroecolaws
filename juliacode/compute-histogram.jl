#= Module for computing the histograms given log-frequency data =#
#/ Start module
module Histogram

#/ Packages
using DataFrames, DataFramesMeta
using FHist 

#################
### FUNCTIONS ###
"""
Compute the histogram

!note: the histogram does not take any cutoff into account, so make sure that the data for
       which the histogram is computed does not contain any data below the cutoff
"""
function compute_fhist(logfrequencies; nbins=16)
    #/ Define binedges
    bmin, bmax = minimum(logfrequencies), maximum(logfrequencies)
    binedges = range(bmin, stop=bmax, length=nbins)
    #/ Compute FHist.Hist1D histogram
    #~ note: put overflow=true to include entries for which logfreq=bmin
    fh = FHist.Hist1D(logfrequencies, binedges=binedges, overflow=true)
    return fh
end

end # module Histogram
#/ End module
