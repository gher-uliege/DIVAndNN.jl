using DIVAnd
using Interpolations
using NCDatasets

include("emodnet_bio_grid.jl");


# on gher17
ds = Dataset(expanduser("~/workspace/ContourExtractor/data/dist2coast.nc"));
Dlon = nomissing(ds["lon"][:])
Dlat = nomissing(ds["lat"][:])

i0,i1 = extrema(findall(gridlon[1]-1 .<= Dlon .<= gridlon[end]+1))
j0,j1 = extrema(findall(gridlat[1]-1 .<= Dlat .<= gridlat[end]+1))

D = nomissing(ds["distance2coast"][i0:i1,j0:j1]);

itp = interpolate((Dlon[i0:i1],Dlat[j1:-1:j0]), D[:,end:-1:1], Gridded(Linear()));
Di = itp(gridlon,gridlat);


datadir = get(ENV,"DATADIR",expanduser("~/tmp/Emodnet-Bio"))


fname = joinpath(datadir,"dist2coast_subset.nc")

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

ncvar = defVar(ds,"distance", Float32, ("lon", "lat")) 
ncvar.attrib["_FillValue"] = Float32(9.96921e36)
ncvar.attrib["missing_value"] = Float32(9.96921e36)
ncvar.attrib["long_name"] = "distance to nearest coast line"


# Define variables

nclon[:] = gridlon
nclat[:] = gridlat
ncvar[:] = Di

close(ds)
