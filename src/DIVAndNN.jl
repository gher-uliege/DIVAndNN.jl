module DIVAndNN

using JSON
using CSV
using DataFrames
using Glob
using Statistics
using Proj4

include("emodnet_bio_loadobs.jl")
include("DIVAnd_covar.jl")
include("emodnet_bio_summary.jl")

export listnames
export loadbyname

end



