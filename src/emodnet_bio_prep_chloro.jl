using Interpolations
using NCDatasets
using PyPlot
using DIVAnd
include("emodnet_bio_grid.jl")

# BalticChla.nc is obtained by subsetting the file complete file
# A20021852018212.L3m_CU_CHL_chlor_a_4km.nc downloaded from
# https://oceancolor.gsfc.nasa.gov/l3/
datafile = joinpath(datadir, "BalticChla.nc")

"""
read_chla_oceancolor(filename, valex)

Read the coordinates and the field from the data file
"""
function read_chla_oceancolor(filename::String, valex::Float64=-999.)
    if isfile(datafile)
        ds = Dataset(filename,"r")
        lon = ds["lon"][:]
        lat = ds["lat"][:][end:-1:1]
        chla = coalesce.(ds["chlor_a"][:][:,end:-1:1], valex)
        chlaf = DIVAnd.ufill(chla, valex);
        close(ds)
    else
        @error "File doesn't exist"
    end

    return lon, lat, chlaf
end


"""
write_chloro_interp(lon, lat, chla_interp, filename)

Write the interpolated field in a new netCDF file
"""
function write_chloro_interp(lon, lat, chla_interp, filename::String)

    Dataset(filename,"c") do ds

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",length(gridlon))
        defDim(ds,"lat",length(gridlat))

        # Define a global attribute
        ds.attrib["title"] = "chlorophyll-a concentration from MODIS-Aqua"

        # Define the variables
        lon = defVar(ds, "lon", Float32, ("lon",))
        lat = defVar(ds, "lat", Float32, ("lat",))
        chla = defVar(ds,"chla", Float32,("lon","lat"))

        lon[:] = gridlon;
        lat[:] = gridlat;
        chla[:,:] = chla_interp;

    end
end

# Read data
@info "Read data"
lon, lat, chla = read_chla_oceancolor(datafile);

# Perform interpolation
@info "Perform interpolation"
itp = interpolate((lon, lat), chla, Gridded(Linear()));
chla_interp = itp(gridlon, gridlat);

# Write in netCDF
@info "Writing new netCDF"
outputfile = joinpath(datadir, "chloro_reinterp.nc");
write_chloro_interp(lon, lat, chla_interp, outputfile);
