#= Module for computing the histograms given log-frequency data =#
#/ Start module
module Histogram

#/ Packages
using FHist
using StatsBase

#################
### FUNCTIONS ###
"""
Compute the histogram

!note: the histogram does not take any cutoff into account, so make sure that the data for
       which the histogram is computed does not contain any data below the cutoff
!note: uses Terrel-Scott rule for the no. of bins, bins = ³√2n, with n the no. of samples
"""
function compute_fhist(logfrequencies; nbins=round(Int, (2*length(logfrequencies))^(1/3)))
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
