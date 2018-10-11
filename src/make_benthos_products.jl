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
function make_analysis(obslon, obslat, g1, g2, g3)

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
    f1ori = exp.(fi1 .+ mean(g1log)) .- 1;
    f2ori = exp.(fi2 .+ mean(g2log)) .- 1;
    f3ori = exp.(fi3 .+ mean(g3log)) .- 1;
    f1ori[f1ori.<0] .= 0.
    f2ori[f2ori.<0] .= 0.
    f3ori[f3ori.<0] .= 0.
    totalfield = f1ori + f2ori + f3ori;

    f1rel, f2rel, f3rel = f1ori./totalfield, f2ori./totalfield, f3ori./totalfield;

    return f1rel, f2rel, f3rel, totalfield;
end;

"""
write_benthos_nc(filename, gridlon, gridlat, field1, field2, field3)

Write the result of the analysis (`make_analysis`) in a netCDF file
"""
function write_benthos_nc(filename::String, gridlon, gridlat,
        field1::Array, field2::Array, field3::Array)

    Dataset(filename,"c") do ds

        nlon = length(gridlon);
        nlat = length(gridlat);

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",nlon);
        defDim(ds,"lat",nlat);

        # Define a global attribute
        ds.attrib["title"] = "Interpolated Benthos "

        # Define the variables and coordinates
        lon = defVar(ds,"lon",Float32,("lon",))
        lat = defVar(ds,"lat",Float32,("lat",))
        g1 = defVar(ds,"g1",Float64,("lon","lat"))
        g2 = defVar(ds,"g2",Float64,("lon","lat"))
        g3 = defVar(ds,"g3",Float64,("lon","lat"))

        # Fill the coord vectors and the fields
        lon[:] = gridlonBenthos;
        lat[:] = gridlatBenthos;

        g1[:,:] = field1,
        g2[:,:] = field2;
        g3[:,:] = field3;

    end
end;

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
fname = joinpath(datadir, "Olivier-Benthos/Benthos_tabSR.csv");
@info "Reading data file $(fname)"
obslon, obslat, g1, g2, g3 = read_benthos(fname);
@info extrema(obslat)
@info extrema(obslon)
@info "Interpolating"
fi1rel, fi2rel, fi3rel, totalfield = make_analysis(obslon, obslat, g1, g2, g3);

outputdir = "../output/"
if !isdir(outputdir)
    @info("Creating output directory $(outputdir)")
    mkdir(outputdir);
end


@info "Write netCDF"
write_benthos_nc(joinpath(outputdir, "Benthos_tabSR.nc"), gridlonBenthos, gridlatBenthos,
fi1rel, fi2rel, fi3rel);

fname = joinpath(datadir, "Olivier-Benthos/Benthos_tabDens.csv");
@info "Reading data file $(fname)"
obslon2, obslat2, g1b, g2b, g3b = read_benthos(fname);
@info extrema(obslat2)
@info extrema(obslon2)
@info "Interpolating"
f1brel, f2brel, f3brel, totalfieldb = make_analysis(obslon2, obslat2, g1b, g2b, g3b);
@info "Write netCDF"
write_benthos_nc(joinpath(outputdir, "Benthos_tabDens.nc"), gridlonBenthos, gridlatBenthos,
fi1rel, fi2rel, fi3rel);
