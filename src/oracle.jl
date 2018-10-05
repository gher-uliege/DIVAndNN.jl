if VERSION >= v"0.7"
    using Test
    using Dates
else
    using Base.Test
end


"""
    read_oracle_file(filename)

Read the coordinates from the Bio-ORACLE layers

Return
    lonr, latr: longitudes and latitudes of the grid ()
    field: the 2-D field storing the variable
    valex: the fill (or exclusion value)

# Example
```julia-repl
julia> lon, lat, temp, valex = read_oracle_file("Present.Surface.Iron.Mean.asc")
```
"""
function read_oracle_file(filename::String)
    if isfile(filename)
        @info "Reading data from file $(filename)"
    else
        @error "File $(filename) does not exist"
        return
    end

    # Open file and read line by line
    open(datafile, "r") do df
        line = readline(df)
        ncols = parse(Int,split(line)[2])
        line = readline(df)
        nrows = parse(Int,split(line)[2])
        line = readline(df)
        XLLCORNER = parse(Float64,split(line)[2])
        line = readline(df)
        YLLCORNER = parse(Float64,split(line)[2])
        line = readline(df)
        CELLSIZE = parse(Float64,split(line)[2])
        line = readline(df)
        valex = parse(Float64, split(line)[2])
        @info("Number of columns (longitude): $(ncols)")
        @info("Number of rows (latitude): $(nrows)")
        @info("Valex = $(valex)")

        # Generate the grid
        lonr = collect(range(XLLCORNER, stop=XLLCORNER + ncols*CELLSIZE, length=ncols));
        latr = collect(range(YLLCORNER, stop=YLLCORNER + nrows*CELLSIZE, length=nrows));

        # Allocate matrix
        field = Array{Float64,2}(undef, nrows, ncols);

        # Read the rest of the data
        iline = 1
        while !eof(df)
            dataline = readline(df)
            field[nrows-iline+1, :] = parse.(Float64,split(dataline));
            iline += 1;
        end

        return lonr, latr, field, valex
    end
end

"""
    subset_oracle(lonr, latr, field, domain)

Subset the coordinates and the fields in the selected domain

Inputs
    lonr, latr: grid coordinates
    field: 2D field storing the variable
    domain: 4-element array with lonmin, lonmax, latmin, latmax

Return
    lond, latd: the grid coordinates in the specified domain
    fieldd: the 2D field in the specified domain

# Example
```julia-repl
julia> lond, latd, temp = subset_oracle(lonr, latr, field., [-5, 30., 35., 50.])
```
"""
function subset_oracle(lonr::Array{Float64,1}, latr::Array{Float64,1},
        field::Array{Float64,2}, domain::typeof(coords))

    goodlon = findall((lonr .>= domain[1]) .& (lonr .<= domain[2]));
    goodlat = findall((latr .>= domain[3]) .& (latr .<= domain[4]));
    lond = lonr[goodlon];
    latd = latr[goodlat];
    fieldd = field[goodlat, goodlon];

    return lond, latd, fieldd
end
