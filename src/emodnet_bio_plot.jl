using PyPlot
using DIVAnd
using OceanPlot
using NCDatasets

if VERSION >= v"0.7"
    using DelimitedFiles
    using Dates
end

include("emodnet_bio_grid.jl")
include("emodnet_bio_loadobs.jl")

outdir = joinpath(datadir,"Results","Zooplankton")
outdir = joinpath("/home/abarth/mnt/nic4/tmp/Emodnet-Bio/","Results","Zooplankton")


fname = joinpath(datadir,"balticZooplankton.csv")
scientificname_accepted = listnames(fname);

sname = scientificname_accepted[1]


#orientation = "vertical"
orientation = "horizontal"


figdir = expanduser("~/Doc/Pres/2018/EMODNET-Bio/Fig3")

mkpath(figdir)

bathname = joinpath(datadir,"gebco_30sec_4.nc");
bathisglobal = true;

blon,blat,bath = DIVAnd.load_bath(bathname,bathisglobal,gridlon,gridlat);
coast() = contourf(blon,blat,bath' .> 0, levels=[0,.5], cmap = "gray")

trans(x) = log(x+1)

function myplot(data,ca)
    pcolor(gridlon,gridlat,copy(trans.(data')));
    clim(ca...)
    colorbar(orientation = orientation);
    coast()
end

XY = DIVAnd.ndgrid(gridlon,gridlat)

@show sname
lon,lat,obstime,value,ids = loadbyname(fname,years,sname)
outname = joinpath(outdir,"DIVAndNN-analysis-$(sname).nc")

ds = Dataset(outname)
value_analysis2 = nomissing(ds["abundance"][:,:,:],NaN)
close(ds)

figure(figsize = (11,7))
for i = 1:length(years)
    syear = years[i]


    sel_year = (Dates.year.(obstime) .== syear)

    xy = (lon[sel_year],lat[sel_year]);
    v = value[sel_year];
    meanv = binobs(xy,v,XY);

    value_analysis2_slice = value_analysis2[:,:,i]

    d = vcat(meanv[isfinite.(meanv)],
             value_analysis2_slice[isfinite.(value_analysis2_slice)])
    ca = extrema(trans.(d))

    clf()
    subplot(1,2,1)
    myplot(value_analysis2_slice,ca)
    title("$(sname) ($(syear))")
    set_aspect_ratio()

    subplot(1,2,2)
    myplot(meanv,ca)
    title("Observations - $(sname) ($(syear))")
    set_aspect_ratio()
    #savefig(joinpath(figdir,"obs-$(sname)-$(syear).png"))
    savefig(joinpath(figdir,"$(sname)-$(syear).png"))
end
