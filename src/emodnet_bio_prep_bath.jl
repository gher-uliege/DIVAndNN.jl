function prep_bath(bathname,bathisglobal,gridlon,gridlat,datadir)

#bathisglobal = true;
#bathname = joinpath(datadir,"gebco_30sec_4.nc");

blon,blat,bath = DIVAnd.load_bath(bathname,bathisglobal,gridlon,gridlat);


@assert blon == gridlon
@assert blat == gridlat

fname = joinpath(datadir,"bathymetry.nc")

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

ncbatymetry = defVar(ds,"batymetry", Float32, ("lon", "lat"))
ncbatymetry.attrib["_FillValue"] = Float32(9.96921e36)
ncbatymetry.attrib["missing_value"] = Float32(9.96921e36)
ncbatymetry.attrib["long_name"] = "bathymetry"


# Define variables

nclon[:] = blon
nclat[:] = blat
ncbatymetry[:] = bath

close(ds)
end
