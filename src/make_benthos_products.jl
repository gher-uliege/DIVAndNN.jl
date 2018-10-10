using DIVAnd
using PyPlot
using NCDatasets
using Missings
using Interpolations
using Plots

if VERSION >= v"0.7"
    using Random
    using DelimitedFiles
    using Statistics
    using Printf
    using FileIO
else
    using Compat: @info, @warn, range, cat
end

include("../src/emodnet_bio_grid.jl");
pyplot()

"""
g1, g2, g3 = read_benthos(filename)

Read the data from the benthos data file
"""
function read_benthos(filename::String)
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

    return obslon, obslat, g1, g2, g3
end


"""
make_analysis(g1, g2, g3)

Perform 1. data transformation,
2. DIVAnd interpolation and
3. inverse transformation
"""
function make_analysis(obslon, obslat, g1::Array, g2::Array, g3::Array)

    # Transformed fields

    g1log = log.(g1.+1);
    g2log = log.(g2.+1);
    g3log = log.(g3.+1);

    # Perform analysis using the selected reference field
    @time fi1,s = DIVAnd.DIVAndrun(mmask[:,:,1],(pm,pn),(xi,yi),
        (obslon,obslat), g1log .- mean(g1log),len,epsilon2,alphabc=2);

    @time fi2,s = DIVAnd.DIVAndrun(mmask[:,:,1],(pm,pn),(xi,yi),
        (obslon,obslat), g2log .- mean(g2log),len,epsilon2,alphabc=2);

    @time fi3,s = DIVAnd.DIVAndrun(mmask[:,:,1],(pm,pn),(xi,yi),
        (obslon,obslat), g3log .- mean(g3log),len,epsilon2,alphabc=2);

    # Tranform back and relative fields
    fi1ori = exp.(fi1 .+ mean(g1log)) .- 1;
    fi2ori = exp.(fi2 .+ mean(g2log)) .- 1;
    fi3ori = exp.(fi3 .+ mean(g3log)) .- 1;
    fi1ori[fi1ori.<0] .= 0.
    fi2ori[fi2ori.<0] .= 0.
    fi3ori[fi3ori.<0] .= 0.
    totalfield = fi1ori + fi2ori + fi3ori;
    fi1rel, fi2rel, fi3rel = fi1ori./totalfield, fi2ori./totalfield, fi3ori./totalfield;

    return f1rel, f2rel, f3real, totalfield
end

# Grid stored in emodnet_bio_grid.jl
xi,yi = DIVAnd.ndgrid(gridlonBenthos, gridlatBenthos);

# Mask
topodir = "/home/ctroupin/Projects/Diva-Workshops/notebooks/data/"
topofile = joinpath(topodir, "gebco_30sec_8.nc");

if isfile(topofile)
    bx, by, b = DIVAnd.load_bath(topofile,true,gridlonBenthos, gridlatBenthos);
    xmask, ymask, mmask = DIVAnd.load_mask(topofile,true,gridlonBenthos, gridlatBenthos,[0]);
else
    @error "Bathymetry file doesn't exist"
end
mmask = mmask[:,:,1]
@info size(mmask)

# Metrics
pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

# Parameter choice: sensitivity analysis
# correlation length
len = 3.;
# obs. error variance normalized by the background error variance
epsilon2 = 2.;


# Read data
fname = joinpath(datadir, "Olivier-Benthos/tab.csv");
@info "Reading data file $(fname)"
obslon, obslat, g1, g2, g3 = read_benthos(fname);
@info extrema(obslat)
@info extrema(obslon)
@info "Interpolating"
f1rel, f2rel, f3rel, totalfield = make_analysis(obslon, obslat, g1, g2, g3);


fname = joinpath(datadir, "Olivier-Benthos/tabsr.csv");
@info "Reading data file $(fname)"
obslon2, obslat2, g1b, g2b, g3b = read_benthos(fname);
@info extrema(obslat2)
@info extrema(obslon2)
@info "Interpolating"
f1brel, f2brel, f3brel, totalfieldb = make_analysis(obslon2, obslat2, g1b, g2b, g3b);
