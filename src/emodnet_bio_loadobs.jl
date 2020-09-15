

struct Format2018
    fname::String
end

"""
    listnames(fname)

List all unique scientific names in the file `fname`
"""
function listnames(df::Format2018)
    fname = df.fname
    data,header = readdlm(fname,',',header = true)
    @assert header == AbstractString["country" "station" "latitude" "longitude" "daycollected" "monthcollected" "yearcollected" "scientificname_accepted" "measurementvalue"]
    header = header[:]

    scientificname_accepted = unique(data[:,findfirst(header .== "scientificname_accepted")])
    return scientificname_accepted
end

"""
    lon,lat,obstime,value,ids = loadbyname(fname,years,scientificname)

Load all data record with the given scientific name and within the year range
"""
function loadbyname(df::Format2018,years,scientificname)
    fname = df.fname

    data,header = readdlm(fname,',',header = true)
    @assert header == AbstractString["country" "station" "latitude" "longitude" "daycollected" "monthcollected" "yearcollected" "scientificname_accepted" "measurementvalue"]
    header = header[:]

    country = Vector{String}(data[:,findfirst(header .== "country")])
    station = Vector{String}(data[:,findfirst(header .== "station")])

    lon = Vector{Float64}(data[:,findfirst(header .== "longitude")])
    lat = Vector{Float64}(data[:,findfirst(header .== "latitude")])

    day = Vector{Int}(data[:,findfirst(header .== "daycollected")])
    month = Vector{Int}(data[:,findfirst(header .== "monthcollected")])
    year = Vector{Int}(data[:,findfirst(header .== "yearcollected")])

    scientificname_accepted = data[:,findfirst(header .== "scientificname_accepted")]

    obstime = DateTime.(year,month,day)

    value = Vector{Float64}(data[:,findfirst(header .== "measurementvalue")])

    sel = (years[1] .<= Dates.year.(obstime)  .<= years[end]) .& (.!ismissing.(value)) .& (scientificname_accepted .== scientificname)
    ids = country .* "-" .* station

    return lon[sel],lat[sel],obstime[sel],value[sel],ids[sel]
end


struct Format2020
    dirname::String
    type::String # analysis or validation
end

function listnames(df::Format2020)
    return sort(unique(map(fn -> rsplit(basename(fn),"-",limit=3)[1],glob("*" * df.type * ".csv",df.dirname))))
end


function loadbyname(df::Format2020,years,scientificname)
    function transform_coords(lon::Array, lat::Array)
        wgs84 = Projection("+proj=longlat +datum=WGS84 +no_defs")
        espgs32361 = Projection("+proj=utm +zone=31 +north +datum=WGS84 +units=m +no_defs")

        npoints = length(lon)
        lonp = Array{Float64, 1}(undef, npoints)
        latp = Array{Float64, 1}(undef, npoints)

        for i = 1:npoints
            lonp[i], latp[i], e = Proj4.transform(espgs32361, wgs84, [lon[i], lat[i], 0.])
        end

        return lonp, latp
    end

    filelist = glob(scientificname * "-*" * df.type * ".csv",df.dirname)
    @assert length(filelist) == 1

    fname = filelist[1]
    df = DataFrame(CSV.read(fname));
    xUTM = df.xUTM
    yUTM = df.yUTM

    lon, lat = transform_coords(xUTM, yUTM);

    obstime = df.date

    value =
        if hasproperty(df, :occurs)
            # Format in April 2020
            df.occurs
        else
            # Format in June 2020
            df.occurrence
        end

    sel = (years[1] .<= Dates.year.(obstime)  .<= years[end]) .& (.!ismissing.(value))
    ids =
        if hasproperty(df, :date_xUTM_yUTM)
            df.date_xUTM_yUTM
        else
            df.eventid
        end

    return lon[sel],lat[sel],obstime[sel],Float64.(value[sel]),ids[sel]
end




function binobs(xy,v,XY)
    sz = size(XY[1])
    n = length(sz)
    sumv = zeros(sz)
    sumv2 = zeros(sz)
    count = zeros(sz)

    xy0 = [XY[i][1] for i = 1:n]

    dxy0 = [XY[i][2,2] - xy0[i] for i = 1:n]

    ind = zeros(Int,n)

    for j in eachindex(v)

        for i = 1:n
            ind[i] = round(Int,(xy[i][j] - xy0[i]) / dxy0[i]) + 1
        end

        count[ind...] += 1
        sumv[ind...] += v[j]
        sumv2[ind...] += v[j]^2
    end

    meanv = sumv./count

    return meanv
end
