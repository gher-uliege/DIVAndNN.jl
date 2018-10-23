using Interpolations
using NCDatasets
using DIVAnd
using Dates
using Statistics
include("emodnet_bio_grid.jl");

function read_emodnet_chem(filename::String,
        varname::String,
        yearmin::Int64, yearmax::Int64, valex::Float64=-999.)

    # Read the data from the seasonal netCDF
    Dataset(filename) do df
        lon = coalesce.(df["lon"][:], NaN);
        lat = coalesce.(df["lat"][:], NaN);
        depth = df["depth"][:];
        time0 = df["time"][:];

        # Select the years of interest
        seldates = findall( (Dates.year.(time0) .≥ yearmin) .&
            (Dates.year.(time0) .≤ yearmax))

        field = df[varname][:];
        @debug size(field)
        field = field[:,:,end, seldates];
        @debug size(field)
        # Replace missing values by exclusion value
        field = coalesce.(field, valex);
        field = DIVAnd.ufill(field, valex);

        return lon, lat, field;
    end
end

"""
list_data_files(datadir)

Provide a list of netCDF files in the directory `datadir`
"""
function list_data_files(datadir::String)::Array
    filelist = [];
    for (root, dirs, files) in walkdir(datadir)
        for file in files
            # List only files ending with nc
            reg = r".nc$"

            if occursin(reg, file)
                fname = joinpath(root, file);
                @debug fname
                push!(filelist, fname)

            end
        end
    end
    return filelist
end

"""
get_mean_field(filelist, years)

Compute the annual mean field based on a list of files and years.
Each file of the list contains the seasonal fields for a series of years.
"""
function get_mean_field(filelist::Array{Any,1}, years::UnitRange{Int64})

    # Read coordinates from 1st file
    df = Dataset(filelist[1])
    nlon = length(df["lon"][:])
    nlat = length(df["lat"][:])
    close(df)

    meanfield = zeros(nlon, nlat, length(years));
    lon = Array{Float64,1}(undef, nlon);
    lat = Array{Float64,1}(undef, nlat);

    for files in filelist
        # Read data from file
        lon, lat, field = read_emodnet_chem(files,
            "Water_body_dissolved_oxygen_concentration",
            minimum(years), maximum(years));

        meanfield = meanfield + field;

    end

    meanfield = meanfield / 4.;

    return lon, lat, meanfield;
end;


"""
write_oxy_interp(gridlon, gridlat, years, field_interp, filename)

Write the longitude, latitude, time (years) and interpolated field in the
file `filename`. 
"""
function write_oxy_interp(gridlon, gridlat, years, field_interp, filename::String)

    Dataset(filename,"c") do ds

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",length(gridlon))
        defDim(ds,"lat",length(gridlat))
        defDim(ds,"time",length(years))

        # Define a global attribute
        ds.attrib["title"] = "Re-interpolated oxygen concentration from EMODnet Chemistry"

        # Define the variables
        lon = defVar(ds, "lon", Float32, ("lon",))
        lat = defVar(ds, "lat", Float32, ("lat",))
        time = defVar(ds, "time", Int64, ("time",))
        oxygen = defVar(ds,"oxygen",Float32,("lon","lat","time"))


        time.attrib["units"] = "Years";
        oxygen.attrib["units"] = "umol/l"

        lon[:] = gridlon;
        lat[:] = gridlat;
        time[:] = collect(years);
        oxygen[:,:] = field_interp;

    end
end

filelist = list_data_files(datadir)
lon, lat, meanfield = get_mean_field(filelist, years);

field_interp = zeros(length(gridlon), length(gridlat), length(years));
xx, yy = DIVAnd.ndgrid(gridlon, gridlat);

# Loop on the selected years
for i = 1:length(years)
    @info "Year: $(years[i])"
    # Perform interpolation
    tmp_itp = interpolate((lon, lat), meanfield[:,:,i], Gridded(Linear()));
    itp =
        @static if VERSION >= v"0.7"
            extrapolate(tmp_itp,Line())
        else
            tmp_itp
        end

    field_interp[:,:,i] = itp.(xx, yy);
end

write_oxy_interp(gridlon, gridlat, years, field_interp, joinpath(datadir, "oxygen_reinterp.nc"));
