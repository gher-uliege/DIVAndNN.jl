import DIVAnd
#using PyPlot
using CSV
using DataFrames
using NCDatasets
using Missings
using Interpolations
using Random
using DelimitedFiles
using Statistics
using Printf
using FileIO
using Dates
using Proj4

include(expanduser("~/src/EMODnet-Biology-Interpolated-Maps/scripts/validate_probability.jl"))
include(expanduser("~/src/EMODnet-Biology-Interpolated-Maps/scripts/PhytoInterp.jl"))

include("DIVAnd_covar.jl")
include("emodnet_bio_loadobs.jl")
include("emodnet_bio_grid.jl")

outdir = joinpath(datadir,"Results","emodnet-bio-2020")

mkpath(outdir)

Random.seed!(1234)

# land-sea mask and domain parameters

maskname = joinpath(datadir,"mask.nc");

ds = Dataset(maskname,"r")
mask = nomissing(ds["mask"][:,:]) .== 1
close(ds)

if ndimensions == 3
    mask = repeat(mask,inner=(1,1,length(years)))
end

bathname = joinpath(datadir,"gebco_30sec_4.nc");
bathisglobal = true;


if ndimensions == 3
    mask2,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat,years);
else
    mask2,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat);
end


# load covariables

covars_fname = [("bathymetry.nc","batymetry",identity),
                ("dist2coast_subset.nc","distance",identity),
                #("Chlorophyll/chloro_reinterp.nc","chla",identity),
                #("oxygen_reinterp2.nc","oxygen",identity),
                ("salinity.nc","salinity",log),
                ("temperature.nc","temperature",identity)
                ]


if ndimensions == 3
    sz = (length(gridlon),length(gridlat),length(years),length(covars_fname)+3)
else
    sz = (length(gridlon),length(gridlat),length(covars_fname)+2)
end

field = zeros(sz)

for i = 1:length(covars_fname)
    (fname,varname,trans) = covars_fname[i]

    Dataset(joinpath(datadir,fname)) do ds
        global tmp
        tmp = nomissing(ds[varname][:],NaN)
        tmp = trans.(tmp)
        if ndimensions == 3
            if ndims(tmp) == 2
                field[:,:,:,i] = repeat(tmp,inner = (1,1,length(years)))
            else
                field[:,:,:,i] = tmp
            end
        else
            field[:,:,i] = tmp
        end
    end
end

X,Y = DIVAnd.ndgrid(gridlon,gridlat)

if ndimensions == 3
    field[:,:,:,end-2] = repeat(X,inner = (1,1,length(years)))
    field[:,:,:,end-1] = repeat(Y,inner = (1,1,length(years)))
    field[:,:,:,end]   = repeat(reshape(years,(1,1,length(years))),inner = (length(gridlon),length(gridlat),1))
else
    field[:,:,end-1] = X
    field[:,:,end] = Y
end

# normalize

for n = 1:size(field,4)
    if ndimensions == 3
        tmp = field[:,:,:,n][mask];
        field[:,:,:,n] = (field[:,:,:,n] .- mean(tmp))./std(tmp)
    else
        tmp = field[:,:,n][mask];
        field[:,:,n] = (field[:,:,n] .- mean(tmp))./std(tmp)
    end
end


#fname = Format2018(joinpath(datadir,"balticZooplankton.csv"))
#fname = Format2020("/home/abarth/src/EMODnet-Biology-Interpolated-Maps/analysis/data/Biddulphia_sinensis1995-2020.csv")

fname = Format2020("/home/abarth/tmp/Emodnet-Bio2020/CSV-split","analysis")

fname_validation = Format2020("/home/abarth/tmp/Emodnet-Bio2020/CSV-split","validation")


scientificname_accepted = listnames(fname);


nameindex = parse(Int,get(ENV,"INDEX","1"))
sname = String(scientificname_accepted[nameindex])
#for sname in scientificname_accepted

@show sname
lon_a,lat_a,obstime_a,value_a,ids_a = loadbyname(fname,years,sname)

lon_cv,lat_cv,obstime_cv,value_cv,ids_cv = loadbyname(fname_validation,years,sname)

#time = Float64.(Dates.year.(obstime_a))


@show value_a[1:min(end,10)]
@show length(value_a)

#=
sum(sel)
plot(lon[sel],lat[sel])
plot(lon[sel],lat[sel],"o")



clf();

scatter(lon[sel],lat[sel],10,value[sel]; cmap = "jet"); colorbar()
=#


#trans,invtrans = DIVAnd.Anam.loglin(Inf, epsilon = 1);
#trans,invtrans = DIVAnd.Anam.notransform();



#=
epsilon2 = 0.5
trans_value = trans.(value)
mvalue = mean(trans_value)
nsamp = 0
varbak,len,dbinfo = DIVAnd.fitlen(
(lon,lat),trans_value,nsamp;
distfun = DIVAnd.distfun_m)

trans_value = trans.(value)
mvalue = mean(trans_value)

f = trans_value .- mvalue

bestfactorl,bestfactore, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D =
DIVAnd.DIVAnd_cv(mask[:,:,1],(pm[:,:,1],pn[:,:,1]),(xi[:,:,1],yi[:,:,1]),(lon,lat),f,len,epsilon2,2,3,0);

#@warn "fix bestfactore"
#bestfactore = 0.2438612845535494

epsilon2 = epsilon2 * bestfactore
=#

len = 50e3
len = 30e3

Random.seed!(1234)

#@show len
#@show epsilon2
#@show bestfactore

Ntries = 1
value_analysis = zeros(size(mask))
value_analysis2 = zeros(size(mask))

lent = 0.6 # years
lent = 0. # years
niter = 100000
#niter = 100000
#niter = 300000
#niter = 10
niter = 2000*100
#niter = 400000 * 16
#testing
#niter = 10
niter = 2000

#for l = 1:Ntries
l=1


#epsilon2ap = 1.5
#epsilon2ap = 2
#epsilon2ap = 5
#epsilon2ap = 50
#@show std(value)
#epsilon2ap = epsilon2
epsilon2ap = 0.5
epsilon2ap = 1.

#NLayers = [size(field)[end],3,1]
NLayers = [size(field)[end],4,1]
#NLayers = [size(field)[end],5,1]
#NLayers = [size(field)[end],2,1]
#NLayers = [size(field)[end],1]
#NLayers = []

learning_rate = 0.00001
learning_rate = 0.001
L2reg = 0.01
dropoutprob = 0.01

xobs_a = if ndimensions == 3
    (lon_a,lat_a,time_a)
else
    (lon_a,lat_a)
end

lenxy = if ndimensions == 3
    (len,len,lent)
else
    (len,len)
end

function plotres(i,lossi,value_analysis2,y)
    vp = validate_probability((gridlon,gridlat),value_analysis2,(lon_cv,lat_cv),value_cv)
	@printf("| %10d | %30.5f | %30.5f | %30.5f |\n",i,lossi,vp,0.)
end

value_analysis2[:],fw0 = analysisprob(
    mask,pmn,xyi,xobs_a,
    value_a,
    lenxy,epsilon2ap,
    field,
    NLayers,
    costfun = nll,
    niter = niter,
    dropoutprob = dropoutprob,
    L2reg = L2reg,
    learning_rate = learning_rate,
	plotres = plotres,
	plotevery = 100
)

#value_analysis2 .= invtrans.(stdf * value_analysis2 .+ mvalue)
#value_analysis2[value_analysis2 .< 0] .= 0

#RMS_diva_covar[l] = validate((gridlon,gridlat,years),value_analysis2,
#		           (lon_cv,lat_cv,time_cv),value_cv)

#return value_analysis,RMS_diva,value_analysis2,RMS_diva_covar

#end

#=

@show mean(RMS_diva), mean(RMS_diva_covar)
@show 1 - mean(RMS_diva_covar)^2/mean(RMS_diva)^2


ncvarattrib = Dict("long_name" => "abundance of $(sname)",
"units" => "1/m2",
"comments" => "number per m2"
)

outname = joinpath(outdir,"DIVAnd-analysis-$(sname).nc")
DIVAnd.save(outname,(gridlon,gridlat,DateTime.(years,1,1)),value_analysis,"abundance"; 
relerr = cpme, ncvarattrib = ncvarattrib)


outname = joinpath(outdir,"DIVAndNN-analysis-$(sname).nc")
DIVAnd.save(outname,(gridlon,gridlat,DateTime.(years,1,1)),value_analysis2,"abundance"; 
relerr = cpme, ncvarattrib = ncvarattrib)

=#


#end

vp = validate_probability((gridlon,gridlat),value_analysis2,(lon_cv,lat_cv),value_cv)
@show vp

outname = joinpath(outdir,"DIVAndNN_$(sname)_interp.nc")

create_nc_results(outname, gridlon, gridlat, value_analysis2, sname;
                  varname = "probability", long_name="occurance probability");

