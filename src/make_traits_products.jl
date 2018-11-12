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
trait_interp = traits_analysis(obslon, obslat, trait)
```
Perform the analysis on `trait` and return the gridded field
"""
function traits_analysis(obslon, obslat, trait)

    # Transformed fields
    traitlog = log.(trait .+ 1);

    # Perform analysis using the selected reference field
    @time fi,s = DIVAnd.DIVAndrun(mmask[:,:,1],(pm,pn),(xi,yi),
        (obslon,obslat), traitlog .- mean(traitlog),len, epsilon2, alphabc=2);

    # Tranform back and relative fields
    trait_interp = exp.(fi .+ mean(traitlog)) .- 1;
    trait_interp[trait_interp.<0] .= 0.

    return trait_interp
end

"""
```julia-repl
trait_interp = traits_error(obslon, obslat, trait)
```
Perform the analysis on `trait` and return the gridded field
"""
function traits_error(obslon, obslat, trait)

    # Transformed fields
    traitlog = log.(trait .+ 1);

    # Perform analysis using the selected reference field
    trait_err = DIVAnd_cpme(mmask[:,:,1],(pm,pn),(xi,yi),(obslon,obslat),
        traitlog .- mean(traitlog),len,epsilon2,alphabc=2);

    return trait_err
end

"""
```julia-repl
make_plot_trait(gridlon, gridlat, field, titletext; vmin=0., vmax=1.)
```
Individual plot for the traits.

"""
function make_plot_trait(gridlon, gridlat, field, titletext; vmin=0., vmax=1.)
    figure()
    ax1 = subplot(1,1,1)
    title(titletext, fontsize=8)
    pcm1 = pcolormesh(gridlon, gridlat, permutedims(field, [2,1]), vmin=vmin, vmax=vmax)
    PyPlot.contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0], colors = [[.5,.5,.5]])
    ax1[:tick_params]("both",labelsize=6)
    gca()[:set_aspect](1/cos(mean([ylim()...]) * pi/180))
    #PyPlot.scatter(obslon, obslat, s=0.5, c=g1log, vmin=0, vmax=10.)
    cb = colorbar(pcm1)[:ax][:tick_params](labelsize=8)
end


"""
```julia-repl
write_traits_nc(filename, gridlon, gridlat, traits3D, traitserr, traitnames)
```
Write the trait analysed and error fields in the netCDF `filename`.

"""
function write_traits_nc(filename::String, gridlon, gridlat,
        traits3D::Array, traitserr::Array, traitnames::Array)

    Dataset(filename,"c") do ds

        nlon = length(gridlon);
        nlat = length(gridlat);
        nt = length(traitnames);

        # Define the dimension "lon" and "lat"
        defDim(ds,"lon",nlon);
        defDim(ds,"lat",nlat);
        defDim(ds,"names",nt)

        # Define a global attribute
        ds.attrib["title"] = "Trait analysis"

        # Define the variables and coordinates
        names = defVar(ds,"names",String,("names",))
        lon = defVar(ds,"lon",Float32,("lon",))
        lat = defVar(ds,"lat",Float32,("lat",))

        # Attributes
        lat.attrib["long_name"] = "Latitude";
        lat.attrib["standard_name"] = "latitude";
        lat.attrib["units"] = "degrees_north";

        lon.attrib["long_name"] = "Longitude";
        lon.attrib["standard_name"] = "longitude";
        lon.attrib["units"] = "degrees_east";

        names.attrib["long_name"] = "Trait name";

        # Fill the coord vectors and the fields
        lon[:] = gridlon;
        lat[:] = gridlat;
        names[:] = traitnames;

        # All the traits are written in the same variable
        # (easier to read). The names of the variables (traits)
        # are stored in the variable
        traits = defVar(ds,"traits",Float64,("lon","lat","names"))
        traits.attrib["long_name"] = "Analysed traits";

        # When using the CPME, the error field is the same for each variable
        traitserror = defVar(ds,"traitserror",Float64,("lon","lat"))
        traitserror.attrib["long_name"] = "Trait error";

        # Fill varibles
        traits[:,:,:] = permutedims(trait_all[:,:,:], [2,3,1]);
        traitserror[:,:] = traitserr;

    end
end;
