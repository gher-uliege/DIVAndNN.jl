using DIVAnd
using Interpolations
using NCDatasets

function prep_mask(bathname,bathisglobal,gridlon,gridlat,years,maskname)

#include("emodnet_bio_grid.jl");

# land-sea mask and domain parameters

mask2,(pm,pn,po),(xi,yi,zi) = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat,years);
blon,blat,bath = DIVAnd.load_bath(bathname,bathisglobal,gridlon,gridlat);
mask = (bath .>= 0);
#mask[1:30,70:end] .= false;
#mask[1:4,1:12] .= false;
#mask[:] = true
X,Y = DIVAnd.ndgrid(gridlon,gridlat)

#=
mask[(X .<= 15.) .& (62 .<= Y)] .= false
mask[(X .<= 10.) .& (Y .<= 54)] .= false
mask[(X .<= 9.3) .& (Y .<= 56)] .= false

mask = DIVAnd.floodfill(mask,CartesianIndex(50,10))
=#

# mask corresponds to the largest area
label = DIVAnd.floodfill(mask)
mask = label .== 1;


ds = Dataset(maskname,"c")
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

ncmask = defVar(ds,"mask", Int8, ("lon", "lat"))
ncmask.attrib["comment"] = "1 sea; 0 land"


# Define variables

nclon[:] = blon
nclat[:] = blat
ncmask[:] = mask

close(ds)

end
