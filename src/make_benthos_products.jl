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
    data,header = readdlm(filename,',',header = true)
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
```julia-repl
f1, f2, f3 = make_analysis(obslon, obslat, g1, g2, g3)
```
Perform
1. data transformation,
2. DIVAnd interpolation and
3. inverse transformation
"""
function make_analysis(obslon, obslat, g1, g2, g3)

    # Transformed fields
    g1log = log.(g1 .+ 1);
    g2log = log.(g2 .+ 1);
    g3log = log.(g3 .+ 1);

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
```julia-repl
err1, err2, err3 = compute_error(obslon, obslat, g1, g2, g3)
```
Compute the error fields using the Clever Poor Man's Estimate.
Data are first transformed using a log function.
"""
function compute_error(obslon, obslat, g1, g2, g3)
    g1log = log.(g1 .+ 1);
    g2log = log.(g2 .+ 1);
    g3log = log.(g3 .+ 1);
    g1_err = DIVAnd_cpme(mmask[:,:,1],(pm,pn),(xi,yi),(obslon,obslat),
        g1log .- mean(g1log),len,epsilon2,alphabc=2);
    g2_err = DIVAnd_cpme(mmask[:,:,1],(pm,pn),(xi,yi),(obslon,obslat),
        g2log .- mean(g2log),len,epsilon2,alphabc=2);
    g3_err = DIVAnd_cpme(mmask[:,:,1],(pm,pn),(xi,yi),(obslon,obslat),
        g3log .- mean(g3log),len,epsilon2,alphabc=2);

    return g1_err, g2_err, g3_err
end

"""
```julia-repl
write_benthos_nc(filename, gridlon, gridlat, field1, field2, field3, err1, err2, err3)
```
Write the result of the analysis (`make_analysis`) and of the error computation (`compute_error`)
in a netCDF file `filename`.
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

        # Interpolated fields
        g1 = defVar(ds,"g1",Float64,("lon","lat"))
        g2 = defVar(ds,"g2",Float64,("lon","lat"))
        g3 = defVar(ds,"g3",Float64,("lon","lat"))

        # Error fields
        g1_err = defVar(ds,"g1_err",Float64,("lon","lat"))
        g2_err = defVar(ds,"g2_err",Float64,("lon","lat"))
        g3_err = defVar(ds,"g3_err",Float64,("lon","lat"))


        # Fill the coord vectors and the fields
        lon[:] = gridlonBenthos;
        lat[:] = gridlatBenthos;

        g1[:,:] = field1;
        g2[:,:] = field2;
        g3[:,:] = field3;

        g1_err[:,:] = err1;
        g2_err[:,:] = err2;
        g3_err[:,:] = err3;

    end
end;

"""
```julia-repl
make_plot_grid(field1, field2, field3, fieldtotal; vmin, vmax)
```

Create a figure with the 3 fields and their total on a 2 by 2 grid.
vmin and vmax are applied as the lower and upper limits for the field values.
"""
function make_plot_grid(field1, field2, field3, fieldtotal; vmin=0, vmax=10.)
    figure()
    ax1 = subplot(2,2,1)
    title("Resistant", fontsize=8)
    pcm1 = pcolormesh(gridlonBenthos, gridlatBenthos, permutedims(field1, [2,1]), vmin=vmin, vmax=vmax)
    PyPlot.contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0],colors = [[.5,.5,.5]])
    ax1[:tick_params]("both",labelsize=6)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    #PyPlot.scatter(obslon, obslat, s=0.5, c=g1log, vmin=0, vmax=10.)
    cb = colorbar(pcm1)[:ax][:tick_params](labelsize=8)
    ax2 = subplot(2,2,2)
    title("Resilient", fontsize=8)
    pcm2 = pcolormesh(gridlonBenthos, gridlatBenthos, permutedims(field2, [2,1]), vmin=vmin, vmax=vmax)
    PyPlot.contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0],colors = [[.5,.5,.5]])
    ax2[:tick_params]("both",labelsize=6)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    #PyPlot.scatter(obslon, obslat, s=0.5, c=g2log, vmin=0, vmax=10.)
    colorbar(pcm2)[:ax][:tick_params](labelsize=8)
    ax3 = subplot(2,2,3)
    title("Vulnerable", fontsize=8)
    pcm3 = pcolormesh(gridlonBenthos, gridlatBenthos, permutedims(field3, [2,1]), vmin=vmin, vmax=vmax)
    PyPlot.contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0],colors = [[.5,.5,.5]])
    ax3[:tick_params]("both",labelsize=6)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    #PyPlot.scatter(obslon, obslat, s=0.5, c=g3log, vmin=0, vmax=10.)
    colorbar(pcm3)[:ax][:tick_params](labelsize=8)
    ax4 = subplot(2,2,4)
    pcm4 = pcolormesh(gridlonBenthos, gridlatBenthos, permutedims(fieldtotal, [2,1]), vmax=300)
    PyPlot.contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0], colors = [[.5,.5,.5]])
    ax4[:tick_params]("both",labelsize=6)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(pcm4)[:ax][:tick_params](labelsize=8)
end
