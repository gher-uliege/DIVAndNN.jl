"""
    listnames(fname)

List all unique scientific names in the file `fname`
"""
function listnames(fname)
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
function loadbyname(fname,years,scientificname)
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



fname = joinpath(datadir,"balticZooplankton.csv")
years = 2007:2013
scientificname = scientificname_accepted[1]


lon,lat,obstime,value,ids = loadbyname(fname,years,scientificname)