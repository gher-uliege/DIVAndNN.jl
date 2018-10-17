using PyPlot
using DIVAnd
using OceanPlot
using NCDatasets
using StatsBase

if VERSION >= v"0.7"
    using DelimitedFiles
    using Dates
end

include("emodnet_bio_grid.jl")
include("emodnet_bio_loadobs.jl")


if get(ENV,"CLUSTER_NAME","") == "nic4"
    outdir = joinpath(datadir,"Results","Zooplankton")
    outdir = joinpath(datadir,"Results","Zooplankton-test4")
else
    outdir = joinpath("/home/abarth/mnt/nic4/tmp/Emodnet-Bio/","Results","Zooplankton")
    outdir = joinpath("/home/abarth/mnt/nic4/tmp/Emodnet-Bio/","Results","Zooplankton-test2")
end

fname = joinpath(datadir,"balticZooplankton.csv")
scientificname_accepted = listnames(fname);


figdir = expanduser("~/Doc/Pres/2018/EMODNET-Bio/Fig3")
figdir = joinpath(outdir,"Fig")

mkpath(figdir)

bathname = joinpath(datadir,"gebco_30sec_4.nc");
bathisglobal = true;

blon,blat,bath = DIVAnd.load_bath(bathname,bathisglobal,gridlon,gridlat);
coast() = contourf(blon,blat,bath' .> 0, levels=[0,.5], cmap = "gray")

trans = x -> log(x+1)
#trans = identity

function myplot(data,ca)
    pcolor(gridlon,gridlat,copy(trans.(data')));
    clim(trans.(ca)...)
    colorbar(orientation = orientation);
    coast()
end

figure(figsize = (7,7))
#orientation = "vertical"
orientation = "horizontal"
XY = DIVAnd.ndgrid(gridlon,gridlat)

sname = scientificname_accepted[1]
for sname in scientificname_accepted

#sname = "Cercopagis (Cercopagis) pengoi"
@show sname
lon,lat,obstime,value,ids = loadbyname(fname,years,sname)
outname = joinpath(outdir,"DIVAndNN-analysis-$(sname).nc")

ds = Dataset(outname)
value_analysis2 = nomissing(ds["abundance"][:,:,:],NaN)
value_analysis2[value_analysis2 .< 0] .= 0
close(ds)

i=1
for i = 1:length(years)
#for i = 1:1
    syear = years[i]


    sel_year = (Dates.year.(obstime) .== syear)

    xy = (lon[sel_year],lat[sel_year]);
    v = value[sel_year];
    meanv = binobs(xy,v,XY);

    value_analysis2_slice = value_analysis2[:,:,i]

    d = vcat(meanv[isfinite.(meanv)],
             value_analysis2_slice[isfinite.(value_analysis2_slice)])

    #d = meanv[isfinite.(meanv)]
    #ca = extrema(d)
    ca = percentile(d[:],[1, 95])

    clf()
#    subplot(1,2,1)
    myplot(value_analysis2_slice,ca)
    title("$(sname) ($(syear))")
    set_aspect_ratio()
    sel = .!isnan.(meanv);
    scatter(XY[1][sel],XY[2][sel],20,trans.(meanv[sel]), edgecolors = "w")
    clim(trans.(ca)...)

    # subplot(1,2,2)
    # myplot(meanv,ca)
    # title("Observations - $(sname) ($(syear))")
    # set_aspect_ratio()
    # #savefig(joinpath(figdir,"obs-$(sname)-$(syear).png"))


    savefig(joinpath(figdir,"$(sname)-$(syear).png"))
end
end
