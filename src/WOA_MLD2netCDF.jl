
using NCDatasets
using DataFrames
using CSV
using PyPlot

"""
read_woa_mld(filename)

Read the content of WOA data file (ascii)
and return the coordinates (lon, lat and mixed layer depth)
"""
function read_woa_mld(filename::String)
    
    lon = Array{Float64,2}(undef,180,360)
    lat = Array{Float64,2}(undef,180,360)
    field = Array{Float64,2}(undef,180,360)
    open(filename) do df
        line = readline(df)
        @info(line)
        iline = 1
        while !eof(df)
            line = readline(df);
            linesplit = split(line);
            lat[iline] = parse(Float64, linesplit[1]);
            lon[iline] = parse(Float64, linesplit[2]);
            field[iline] = parse(Float64, linesplit[3]);
            iline += 1
        end    
    end
    lon2write = sort(unique(lon));
    lat2write = sort(unique(lat));
    return lon2write, lat2write, field
end

"""
write_MLD_netCDF(outputfile, datadir)

Create a netCDF file storing the MLD from the World Ocean Atlas

Inputs
    outputfile: name of the new netCDF file
    datadir: directory where the WOA ascii files (mld.pt.0mm.001) are located
"""

function write_MLD_netCDF(outputfile::String, datadir="./MLD/WOA/"::String)

    nlon = 360
    nlat = 180
    nt = 12
    valex = -99.999

    Dataset(outputfile, "c") do ds

        # Dimensions
        defDim(ds,"lon",nlon)
        defDim(ds,"lat",nlat)
        defDim(ds,"time",nt)

        # Define a global attribute
        ds.attrib["title"] = "MLD monthly Climatology from World Ocean Atlas"
        ds.attrib["source"] = "https://www.nodc.noaa.gov/OC5/WOA94/mix.html"

        # Variables
        lon = defVar(ds,"lon",Float64,("lon",))
        lat = defVar(ds,"lat",Float64,("lat",))
        time0 = defVar(ds,"time",Float64,("time",))

        lat.attrib["standard_name"] = "latitude"
        lat.attrib["long_name"] = "Latitude coordinate"
        lat.attrib["units"] = "degrees_north"
        lat.attrib["valid_min"] = -90.
        lat.attrib["valid_max"] = 90. 
        lat.attrib["axis"] = "Y"

        lon.attrib["standard_name"] = "longitude"
        lon.attrib["long_name"] = "Longitude coordinate"
        lon.attrib["units"] = "degrees_east"
        lon.attrib["valid_min"] = -180.
        lon.attrib["valid_max"] = 180. 
        lon.attrib["axis"] = "X"

        MLD = defVar(ds,"MLD",Float64,("lon","lat","time"))
        MLD.attrib["units"] = "meters"
        MLD.attrib["standard_name"] = "ocean_mixed_layer_thickness_defined_by_temperature"
        MLD.attrib["coordinates"] = "time lat lon";
        MLD.attrib["missing_value"] = valex;

        for i=1:12
            mmm = lpad(i,3,"0")
            datafile = joinpath(datadir, "mld.pt.$(mmm).001")
            if isfile(datafile)
                @info("Reading file $(datafile)")
            else
                @warn("File $(datafile) doesn't exist")
            end

            lon2write, lat2write, field = read_woa_mld(datafile);

            time0[i] = i
            if i == 1
                lon[:] = lon2write;
                lat[:] = lat2write;
            end
            # Modify field to have lon between -180. to 180.
            field = permutedims(field, [2,1]);
            MLD[:,:,i] = vcat(field[181:360,:],field[1:180,:]);
        end
    end
end

write_MLD_netCDF("pp.nc", "../Datasets/MLD/WOA/")


