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

"""
```julia-repl
obslon, obslat, obsyear, g1, g2, g3, g4 = read_fish(filename)
```
Read the coordinates, the year and the abundance of each of the 4 types of fish
"""
function read_fish(filename::String)
    data,header = readdlm(filename,',',header = true)
    header = header[:]
    if "year" in header
        @info "Working on a temporal data file"
        stationname = Vector{String}(data[:,findfirst(header .== "samp")]);
        obsyear = Vector{Int32}(data[:,findfirst(header .== "year")]);
    else
        @info "Working on a spatial data file"
        obsyear = undef;
    end;

    obslon = Vector{Float64}(data[:,findfirst(header .== "x")]);
    obslat = Vector{Float64}(data[:,findfirst(header .== "y")]);

    g1 = Vector{Float64}(data[:,findfirst(header .== "g1")]);
    g2 = Vector{Float64}(data[:,findfirst(header .== "g2")]);
    g3 = Vector{Float64}(data[:,findfirst(header .== "g3")]);
    g4 = Vector{Float64}(data[:,findfirst(header .== "g4")]);

    @info "Number of data points: $(length(g1))"

    return obslon, obslat, obsyear, g1, g2, g3, g4
end

"""
```julia-repl
add_mask(bx, by, b)
```
Add the land-sea mask to the plot with a grey (0.5) color
"""
function add_mask(bx, by, b)
    PyPlot.contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0], colors = [[.5,.5,.5]])
end

"""
```julia-repl
plot_fish_data(obslonS, obslatS, g1, g2, g3, g4)
```
Make a scatter of the data values, one subplot per fish category
"""
function plot_fish_data(obslonS, obslatS, g1, g2, g3, g4)
    figure()
    ax1 = subplot(2,2,1)
    ax1[:tick_params]("both",labelsize=6)
    scat1 = PyPlot.scatter(obslonS, obslatS, s=.01, c=g1)
    add_mask(bx, by, b)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(scat1)[:ax][:tick_params](labelsize=8)
    ax2 = subplot(2,2,2)
    ax2[:tick_params]("both",labelsize=6)
    scat2 = PyPlot.scatter(obslonS, obslatS, s=.01, c=g2)
    add_mask(bx, by, b)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(scat2)[:ax][:tick_params](labelsize=8)
    ax3 = subplot(2,2,3)
    ax3[:tick_params]("both",labelsize=6)
    scat3 = PyPlot.scatter(obslonS, obslatS, s=.01, c=g3)
    add_mask(bx, by, b)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(scat3)[:ax][:tick_params](labelsize=8)
    ax4 = subplot(2,2,4)
    ax4[:tick_params]("both",labelsize=6)
    scat4 = PyPlot.scatter(obslonS, obslatS, s=.01, c=g4)
    add_mask(bx, by, b)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(scat4)[:ax][:tick_params](labelsize=8)
end

"""
```julia-repl
write_fish_nc(filename, gridlon, gridlat, field1, field2, field3, field4,
             err1, err2, err3, err4)
```
Write the interpolated and error fields for the 4 fish categories in the
netCDF file `filename`
"""
function write_fish_nc(filename::String, gridlon, gridlat,
        field1::Array, field2::Array, field3::Array, field4::Array,
        err1::Array, err2::Array, err3::Array, err4::Array)

    Dataset(filename,"c") do ds

        nlon = length(gridlon);
        nlat = length(gridlat);

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",nlon);
        defDim(ds,"lat",nlat);

        # Define a global attribute
        ds.attrib["title"] = "Interpolated Fish"

        # Define the variables and coordinates
        lon = defVar(ds,"lon",Float32,("lon",))
        lat = defVar(ds,"lat",Float32,("lat",))

        # Attributes
        lat.attrib["long_name"] = "Latitude";
        lat.attrib["standard_name"] = "latitude";
        lat.attrib["units"] = "degrees_north";

        lon.attrib["long_name"] = "Longitude";
        lon.attrib["standard_name"] = "longitude";
        lon.attrib["units"] = "degrees_east";

        # Interpolated fields
        g1 = defVar(ds,"g1",Float64,("lon","lat"))
        g2 = defVar(ds,"g2",Float64,("lon","lat"))
        g3 = defVar(ds,"g3",Float64,("lon","lat"))
        g4 = defVar(ds,"g4",Float64,("lon","lat"))

        # Error fields
        g1_err = defVar(ds,"g1_err",Float64,("lon","lat"))
        g2_err = defVar(ds,"g2_err",Float64,("lon","lat"))
        g3_err = defVar(ds,"g3_err",Float64,("lon","lat"))
        g4_err = defVar(ds,"g4_err",Float64,("lon","lat"))

        # Fill the coord vectors and the fields
        lon[:] = gridlon;
        lat[:] = gridlat;

        g1[:,:] = field1;
        g2[:,:] = field2;
        g3[:,:] = field3;
        g4[:,:] = field4;

        g1_err[:,:] = err1;
        g2_err[:,:] = err2;
        g3_err[:,:] = err3;
        g4_err[:,:] = err4;

    end

end;

"""
```julia-repl
plot_fish_results(gridlon, gridlat, f1, f2, f3, f4, bx, by, b)
```
Make the plot of the 4 interpolated fields f1, f2, f3, f4
and overlay the mask based on the bathymetry defined by bx, by and b
"""
function plot_fish_results(gridlon, gridlat, f1, f2, f3, f4, bx, by, b)
    figure()
    ax1 = subplot(2,2,1)
    ax1[:tick_params]("both",labelsize=6)
    pcm1 = PyPlot.pcolormesh(gridlon, gridlat, permutedims(f1, [2,1]), vmin=0, vmax=1.)
    add_mask(bx, by, b)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(pcm1)[:ax][:tick_params](labelsize=8)
    ax2 = subplot(2,2,2)
    ax2[:tick_params]("both",labelsize=6)
    pcm2 = PyPlot.pcolormesh(gridlon, gridlat, permutedims(f2, [2,1]), vmin=0, vmax=1.)
    add_mask(bx, by, b)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(pcm2)[:ax][:tick_params](labelsize=8)
    ax3 = subplot(2,2,3)
    ax3[:tick_params]("both",labelsize=6)
    pcm3 = PyPlot.pcolormesh(gridlon, gridlat, permutedims(f3, [2,1]), vmin=0, vmax=1.)
    add_mask(bx, by, b)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(pcm3)[:ax][:tick_params](labelsize=8)
    ax4 = subplot(2,2,4)
    ax4[:tick_params]("both",labelsize=6)
    pcm4 = PyPlot.pcolormesh(gridlon, gridlat, permutedims(f4, [2,1]), vmin=0, vmax=1.)
    add_mask(bx, by, b)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    colorbar(pcm4)[:ax][:tick_params](labelsize=8)
end

"""
```julia-repl
write_fish_time_nc(filename, gridlon, gridlat, gridtime,
                    field1, field2, field3, field4,
                    err1, err2, err3, err4)
```
Write a netCDF file containing the gridded and error fields for the
temporal fish product
"""
function write_fish_time_nc(filename::String, gridlon, gridlat, gridtime,
        field1::Array, field2::Array, field3::Array, field4::Array,
        err1::Array, err2::Array, err3::Array, err4::Array)

    Dataset(filename,"c") do ds

        nlon = length(gridlon);
        nlat = length(gridlat);
        ntimes = length(gridtime);

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",nlon);
        defDim(ds,"lat",nlat);
        defDim(ds,"time",ntimes);

        # Define a global attribute
        ds.attrib["title"] = "Interpolated Fish - temporal product"

        # Define the variables and coordinates
        lon = defVar(ds,"lon",Float32,("lon",))
        lat = defVar(ds,"lat",Float32,("lat",))
        time = defVar(ds,"time",Float32,("time",))

        # Attributes
        lat.attrib["long_name"] = "Latitude";
        lat.attrib["standard_name"] = "latitude";
        lat.attrib["units"] = "degrees_north";

        lon.attrib["long_name"] = "Longitude";
        lon.attrib["standard_name"] = "longitude";
        lon.attrib["units"] = "degrees_east";

        time.attrib["long_name"] = "Time";
        time.attrib["standard_name"] = "time";
        time.attrib["units"] = "years";

        # Interpolated fields
        g1 = defVar(ds,"g1",Float64,("lon","lat","time"))
        g2 = defVar(ds,"g2",Float64,("lon","lat","time"))
        g3 = defVar(ds,"g3",Float64,("lon","lat","time"))
        g4 = defVar(ds,"g4",Float64,("lon","lat","time"))

        # Error fields
        g1_err = defVar(ds,"g1_err",Float64,("lon","lat","time"))
        g2_err = defVar(ds,"g2_err",Float64,("lon","lat","time"))
        g3_err = defVar(ds,"g3_err",Float64,("lon","lat","time"))
        g4_err = defVar(ds,"g4_err",Float64,("lon","lat","time"))

        # Fill the coord vectors and the fields
        lon[:] = gridlon;
        lat[:] = gridlat;
        time[:] = gridtime;

        g1[:,:,:] = field1;
        g2[:,:,:] = field2;
        g3[:,:,:] = field3;
        g4[:,:,:] = field4;

        g1_err[:,:,:] = err1;
        g2_err[:,:,:] = err2;
        g3_err[:,:,:] = err3;
        g4_err[:,:,:] = err4;
    end

end;

"""
```julia-repl
plot_fish_results_time(years,gridlon, gridlat, f1, bx, by, b)
```
Create the general plot for the fish products.
One plot per fish type, all the years on the plot.
"""
function plot_fish_results_time(years,gridlon, gridlat, f1, bx, by, b)

    yearinit = collect(years) .-1;
    yearinit[1] = yearinit[2];
    yearend = collect(years) .+ 1;
    yearend[end] = yearend[end-1]

    ffig = figure("fish_time",figsize=(12,12))
    for ii = 1:size(f1)[1]
        ax = subplot(5,4,ii)
        title("$(yearinit[ii])-$(yearend[ii])", fontsize=8)
        ax[:tick_params]("both",labelsize=6)
        pcm = PyPlot.pcolormesh(gridlon, gridlat, permutedims(f1[ii,:,:], [2,1]), vmin=0, vmax=1.)
        add_mask(bx, by, b)
        gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
        colorbar(pcm)[:ax][:tick_params](labelsize=8)
    end
end
