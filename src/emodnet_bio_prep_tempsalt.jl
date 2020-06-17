import DIVAnd
#using PyPlot
#using CSV
using NCDatasets
using Missings
using Interpolations
using Base.Threads

include("emodnet_bio_grid.jl");


for (fname,varname) in data_TS
    @show fname
    ds = Dataset(fname)
    k_index = 1
    S = ds[varname][:,:,k_index,:];
    lon = nomissing(ds["lon"][:])
    lat = nomissing(ds["lat"][:])
    SS = nomissing(S,NaN);
    close(ds)

    @info "skip time instance without data"
    mask = .!isnan.(SS);
    count = sum(sum(mask,dims=1),dims=2)[:]
    n = findall(count .> 0)
    SS = SS[:,:,n]
    mask = mask[:,:,n]

    @info "fill"
    SS = DIVAnd.ufill(SS,mask)
#=    Threads.@threads for k = 1:size(mask,3)
        @show k
        SS[:,:,k] = DIVAnd.ufill(SS[:,:,k],mask[:,:,k])
    end
=#
    @info "average"
    SS2 = mean(SS,dims = 3)[:,:,1]


    @info "interpolate"
    itp = interpolate((lon,lat), SS2, Gridded(Linear()));
    SSi = itp(gridlon,gridlat);

    @show extrema(SSi)
    fname = joinpath(datadir,"$(lowercase(varname)).nc")

    ds = Dataset(fname,"c")
    # Dimensions

    ds.dim["lon"] = length(gridlon)
    ds.dim["lat"] = length(gridlat)

    # Declare variables

    nclon = defVar(ds,"lon", Float64, ("lon",))
    nclon.attrib["units"] = "degrees_east"
    nclon.attrib["standard_name"] = "longitude"
    nclon.attrib["long_name"] = "longitude"

    nclat = defVar(ds,"lat", Float64, ("lat",))
    nclat.attrib["units"] = "degrees_north"
    nclat.attrib["standard_name"] = "latitude"
    nclat.attrib["long_name"] = "latitude"

    ncvar = defVar(ds,lowercase(varname), Float32, ("lon", "lat"))
    ncvar.attrib["_FillValue"] = Float32(9.96921e36)
    ncvar.attrib["missing_value"] = Float32(9.96921e36)
    ncvar.attrib["long_name"] = varname


    # Define variables

    nclon[:] = gridlon
    nclat[:] = gridlat
    ncvar[:] = SSi

    close(ds)


end
