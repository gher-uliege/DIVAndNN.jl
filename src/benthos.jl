
include("emodnet_bio_loadobs.jl")


fname = joinpath(datadir,"Olivier-Benthos/tab.csv");
data,header = readdlm(fname,',',header = true)
header = header[:]

# "data","x","y","sta","g1","g2","g3"
dataname = Vector{String}(data[:,findfirst(header .== "data")]);

obslon = Vector{Float64}(data[:,findfirst(header .== "x")]);
obslat = Vector{Float64}(data[:,findfirst(header .== "y")]);

stationname = Vector{String}(data[:,findfirst(header .== "sta")]);

g1 = Vector{Float64}(data[:,findfirst(header .== "g1")]);
g2 = Vector{Float64}(data[:,findfirst(header .== "g2")]);
g3 = Vector{Float64}(data[:,findfirst(header .== "g3")]);
