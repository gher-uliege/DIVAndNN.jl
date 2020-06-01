import DIVAnd
#using PyPlot
#using CSV
using NCDatasets
using Missings
using Interpolations
using Random
using DelimitedFiles
using Statistics
using Printf
using FileIO
using Dates


function interpcv(xyi,value_analysis,xy)
    ncv = length(xy[1])
    value_analysis_cv = zeros(ncv)

    xyi = map(x -> collect(x),xyi);

    tmp_itp = interpolate(xyi, value_analysis, Gridded(Linear()))
    itp = extrapolate(tmp_itp,Line())

    #save("tmp.jld2","xyi",xyi,"value_analysis",value_analysis,"xy",xy)
    xyi = zeros(length(xy))

    for i = 1:ncv
    	for j = 1:length(xy)
	        xyi[j] = xy[j][i]
        end

        value_analysis_cv[i] = itp(xyi...)
    end

    return value_analysis_cv
end

function validate(xyi,value_analysis,xy_cv,value_cv )
    if isempty(value_cv)
        return NaN
    end

    #@show extrema(value_analysis[.!isnan.(value_analysis)])
    value_analysis_cv = interpcv(xyi,value_analysis,xy_cv)

    sels = .!isnan.(value_analysis_cv);
    #@show std(value_cv[sels] - value_analysis_cv[sels])
    RMS = sqrt(mean(abs2,value_cv[sels] - value_analysis_cv[sels]))
    #@show RMS
    return RMS
end


include("DIVAnd_covar.jl")
include("emodnet_bio_grid.jl")
include("emodnet_bio_loadobs.jl")

outdir = joinpath(datadir,"Results","Zooplankton")
outdir = joinpath(datadir,"Results","Zooplankton-test5")
mkpath(outdir)

Random.seed!(1234)


# land-sea mask and domain parameters

maskname = joinpath(datadir,"mask.nc");

ds = Dataset(maskname,"r")
mask = nomissing(ds["mask"][:,:]) .== 1
close(ds)
mask = repeat(mask,inner=(1,1,length(years)))
bathname = joinpath(datadir,"gebco_30sec_4.nc");
bathisglobal = true;
mask2,(pm,pn,po),(xi,yi,zi) = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat,years);


# load covariables

covars_fname = [("bathymetry.nc","batymetry",identity),
                ("dist2coast_subset.nc","distance",identity),
                ("Chlorophyll/chloro_reinterp.nc","chla",identity),
                ("oxygen_reinterp2.nc","oxygen",identity),
                ("salinity.nc","salinity",log),
                ("temperature.nc","temperature",identity)]

sz = (length(gridlon),length(gridlat),length(years),length(covars_fname)+3)
field = zeros(sz)

for i = 1:length(covars_fname)
    (fname,varname,trans) = covars_fname[i]

    Dataset(joinpath(datadir,fname)) do ds
        global tmp
        tmp = nomissing(ds[varname][:],NaN)
        tmp = trans.(tmp)
        if ndims(tmp) == 2
            field[:,:,:,i] = repeat(tmp,inner = (1,1,length(years)))
        else
            field[:,:,:,i] = tmp
        end
    end
end

X,Y = DIVAnd.ndgrid(gridlon,gridlat)

field[:,:,:,end-2] = repeat(X,inner = (1,1,length(years)))
field[:,:,:,end-1] = repeat(Y,inner = (1,1,length(years)))
field[:,:,:,end]   = repeat(reshape(years,(1,1,length(years))),inner = (length(gridlon),length(gridlat),1))

# normalize

for n = 1:size(field,4)
    tmp = field[:,:,:,n][mask];
    field[:,:,:,n] = (field[:,:,:,n] .- mean(tmp))./std(tmp)
end


fname = joinpath(datadir,"balticZooplankton.csv")
scientificname_accepted = listnames(fname);


nameindex = parse(Int,get(ENV,"INDEX","1"))
sname = scientificname_accepted[nameindex]
#for sname in scientificname_accepted

    @show sname
    lon,lat,obstime,value,ids = loadbyname(fname,years,sname)
 
    time = Float64.(Dates.year.(obstime))

    @show value[1:min(end,10)]
    @show length(value)

#=
sum(sel)
plot(lon[sel],lat[sel])
plot(lon[sel],lat[sel],"o")



clf();

scatter(lon[sel],lat[sel],10,value[sel]; cmap = "jet"); colorbar()
=#





trans,invtrans = DIVAnd.Anam.loglin(Inf, epsilon = 1);
#trans,invtrans = DIVAnd.Anam.notransform();

if false
    pcolor(field[:,:,1]')



    PyPlot.plt[:hist](value[sel],100)
    title("Histogram");
    savefig(expanduser("~/tmp/Emodnet-Bio/hist_value.png"))



    clf();
    PyPlot.plt[:hist](trans.(value[sel]),100)
    title("Histogram (transformed)");
    savefig(expanduser("~/tmp/Emodnet-Bio/hist_value_trans.png"))

end


#=
XY = DIVAnd.ndgrid(gridlon,gridlat)

xy = (lon[sel],lat[sel],time[sel]);
v = value[sel];
meanv = binobs(xy,v,XY);


=#
Ntries = 1


RMS_diva = zeros(Ntries)
RMS_diva_covar = zeros(Ntries)



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

Random.seed!(1234)

@show len
@show epsilon2
@show bestfactore

Ntries = 1
value_analysis = zeros(size(mask))
value_analysis2 = zeros(size(mask))

lent = 0.6 # years
lent = 0. # years
niter = 100000
#niter = 100000
#niter = 300000
#niter = 10
#niter = 2000
niter = 400000 * 16
#testing
niter = 10

#for l = 1:Ntries
    l=1

    # 5% for cross-validation
    selcv = rand(Float64,size(value)) .< 0.10
    @show sum(selcv)

    #    value_analysis,RMS_diva[l],value_analysis2,RMS_diva_covar[l] = comp(selcv)

    lon_a = lon[.!selcv]
    lat_a = lat[.!selcv]
    time_a = time[.!selcv]
    value_a = Vector{Float64}(value[.!selcv])

    lon_cv = lon[selcv]
    lat_cv = lat[selcv]
    time_cv = time[selcv]
    value_cv = Vector{Float64}(value[selcv])

    trans_value_a = trans.(value_a)
    mvalue = mean(trans_value_a)

    f = trans_value_a .- mvalue

    value_analysis[:],s = DIVAnd.DIVAndrun(
        mask,(pm,pn,po),(xi,yi,zi),(lon_a,lat_a,time_a),
        f,
        (len,len,lent),epsilon2;
        alphabc = 0,
    );

    cpme = DIVAnd.DIVAnd_cpme(mask,(pm,pn,po),(xi,yi,zi),(lon_a,lat_a,time_a),f,(len,len,lent),epsilon2)

    value_analysis .= invtrans.(value_analysis .+ mvalue);

    #lon,lat,10,value
    RMS_diva[l] = validate((gridlon,gridlat,years),value_analysis,
    		           (lon_cv,lat_cv,time_cv),value_cv)

    #epsilon2ap = 0.4
    #epsilon2ap = 1.5
    epsilon2ap = 2
    #epsilon2ap = 5
    #epsilon2ap = 50
    #@show std(value)
    epsilon2ap = epsilon2

    #NLayers = [size(field)[end],3,1]
    NLayers = [size(field)[end],4,1]
    #NLayers = [size(field)[end],5,1]
    #NLayers = [size(field)[end],2,1]
    #NLayers = [size(field)[end],1]
    #NLayers = []
    RMS_diva_covar_his = Float64[]
    stdf = std(f)

    value_analysis2[:],fw0 = analysisprob(
        mask,(pm,pn,po),(xi,yi,zi),(lon_a,lat_a,time_a),
        f / stdf,
        (len,len,lent),epsilon2ap,
        field,
        NLayers,
        costfun = regression,
        niter = niter,
        dropoutprob = 0.01,
        L2reg = 0.01,
        learning_rate = 0.0001,
	    plotres = (i,lossi,field,y) -> begin
            #@show i
	        value_analysis2 = invtrans.(stdf * field .+ mvalue)
                value_analysis2[value_analysis2 .< 0] .= 0
    	    RMS = validate((gridlon,gridlat,years),value_analysis2,
    	    	           (lon_cv,lat_cv,time_cv),value_cv)
            push!(RMS_diva_covar_his,RMS)
	        @printf("| %10d | %30.5f | %30.5f | %30.5f |\n",i,lossi,RMS,1 - RMS^2/mean(RMS_diva)^2)
	    end,
	    plotevery = 1000
    )

    value_analysis2 .= invtrans.(stdf * value_analysis2 .+ mvalue)
    value_analysis2[value_analysis2 .< 0] .= 0

    RMS_diva_covar[l] = validate((gridlon,gridlat,years),value_analysis2,
    		           (lon_cv,lat_cv,time_cv),value_cv)

    #return value_analysis,RMS_diva,value_analysis2,RMS_diva_covar

#end



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




#end



