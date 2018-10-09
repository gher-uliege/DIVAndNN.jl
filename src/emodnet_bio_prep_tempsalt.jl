import DIVAnd
#using PyPlot
#using CSV
using NCDatasets
using Missings
using Interpolations


include("emodnet_bio_grid.jl");


global SSi

for (fname,varname) in [("http://sdn.oceanbrowser.net:8081/data/SeaDataNet-domains/Baltic/Salinity.19002012.4Danl.nc","Salinity")
                        ("http://sdn.oceanbrowser.net:8081/data/SeaDataNet-domains/Baltic/Temperature.19002012.4Danl.nc","Temperature")]


    ds = Dataset(fname)
    S = ds[varname][:,:,end,:];
    lon = nomissing(ds["lon"][:])
    lat = nomissing(ds["lat"][:])
    SS = nomissing(S,NaN);
    close(ds)

    
    SS = DIVAnd.ufill(SS,.!isnan.(SS))
    
    SS2 = mean(SS,3)[:,:,1]
    
    
    itp = interpolate((lon,lat), SS2, Gridded(Linear()));
    SSi = itp(gridlon,gridlat);
    


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
