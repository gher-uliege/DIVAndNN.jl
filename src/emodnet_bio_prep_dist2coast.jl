using DIVAnd
using Interpolations
using NCDatasets

#=
include("emodnet_bio_grid.jl");
# on gher17
fname_dist2coast = expanduser("~/workspace/ContourExtractor/data/dist2coast.nc")

datadir = get(ENV,"DATADIR",expanduser("~/tmp/Emodnet-Bio"))
interp_fname = joinpath(datadir,"dist2coast_subset.nc")

prep_dist2coast(fname_dist2coast,gridlon,gridlat,interp_fname; extra = 1, varname = "distance2coast")
=#


"""
    DIVAndNN.prep_dist2coast(fname,gridlon,gridlat,interp_fname; extra = 1, varname = "dist")

Interpolates the dist2coast_1deg dataset (https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/dist2coast_1deg)
onto the grid defined by `(gridlon,gridlat)` and saves the result in the file `interp_fname`.

## Example

```julia
gridlon = -90.:0.5:40.
gridlat = 30.:0.5:80.

interp_fname = joinpath(datadir,"dist2coast_subset.nc")
fname_dist2coast = "https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/dist2coast_1deg"

DIVAndNN.prep_dist2coast(fname_dist2coast,gridlon,gridlat,interp_fname)
```
"""
function prep_dist2coast(fname_dist2coast,gridlon,gridlat,interp_fname; extra = 1, varname = "dist")
    ds = Dataset(fname_dist2coast);
    Dlon = ds["lon"][:]
    Dlat = ds["lat"][:]

    i0,i1 = extrema(findall(gridlon[1]-extra .<= Dlon .<= gridlon[end]+extra))
    j0,j1 = extrema(findall(gridlat[1]-extra .<= Dlat .<= gridlat[end]+extra))

    if fname_dist2coast == "https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/dist2coast_1deg"
        # special case 
        # _FillValue is zero and not missing
        D = ds[varname].var[i0:i1,j0:j1]
    else
        D = nomissing(ds[varname][i0:i1,j0:j1])
    end

    close(ds)

    lon = Dlon[i0:i1]
    lat = Dlat[j0:j1]

    if lat[2] < lat[1]
        lat = reverse(lat)
        D = reverse(D,dims=2)
    end

    varnamei = "distance"
    DIVAndNN.saveinterp((lon,lat),D,(gridlon,gridlat),varnamei,interp_fname)
end
