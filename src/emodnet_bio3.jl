import DIVAnd
#using PyPlot
#using CSV
using NCDatasets
using Missings
using Interpolations
using MAT

if VERSION >= v"0.7"
    using Random
    using DelimitedFiles
    using Statistics
    using Printf
    using FileIO
else
    using Compat: @info, @warn, range, cat
end

function mad(x)
    sel = .!isnan.(x)
    return median(abs.(x[sel] - median(x[sel])))
end

function interpcv(xyi,value_analysis,xy)
    value_analysis_cv = zeros(size(lon_cv))

    @static if VERSION >= v"0.7"
        # https://github.com/JuliaMath/Interpolations.jl/issues/248
        xyi = map(x -> collect(x),xyi);
    end

    tmp_itp = interpolate(xyi, value_analysis, Gridded(Linear()))

    itp =
        @static if VERSION >= v"0.7"
            extrapolate(tmp_itp,Line())
        else
            tmp_itp
        end

    #save("tmp.jld2","xyi",xyi,"value_analysis",value_analysis,"xy",xy)
    xyi = zeros(length(xy))

    for i = 1:length(lon_cv)
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


function surfacemean(fname,varname)
    step = 2

    ds = Dataset(fname)
    S = ds[varname][:,:,end,:];
    lon = nomissing(ds["lon"][1:step:end])
    lat = nomissing(ds["lat"][1:step:end])
    SS = nomissing(S);
    SS[ismissing.(S)] .= NaN;
    close(ds)
    return lon,lat,mean(SS,3)[1:2:end,1:2:end,1]
end

coast() = contourf(bx,by,b' .> 0, levels=[0,.5], cmap = "gray")

include("DIVAnd_covar.jl")

if VERSION >= v"0.7"
    Random.seed!(1234)
else
    srand(1234)
end

datadir = get(ENV,"DATADIR",expanduser("~/tmp/Emodnet-Bio"))

bathname = joinpath(datadir,"gebco_30sec_4.nc");


bathisglobal = true;

#df = CSV.read(joinpath(datadir,"20180410_232725_15acd2c3d3150f.csv"))
# measurement
#df = CSV.read(joinpath(datadir,"20180503_144507_15aeb04539abf4.csv"))

data,header = readdlm(joinpath(datadir,"outalex.csv"),',',header = true)
header = header[:]

lon = Vector{Float64}(data[:,findfirst(header .== "decimalLongitude")])
lat = Vector{Float64}(data[:,findfirst(header .== "decimalLatitude")])
sname = "Acartia bifilosa"


value = Vector{Float64}(data[:,findfirst(header .== sname)])
obstime = DateTime.(data[:,findfirst(header .== "date")])
time = Float64.(Dates.year.(obstime))


#year = 2007;
#year = 2010;
#years = 2000:2012
years = 2007:2013

#for year in years

sel = (years[1] .<= Dates.year.(obstime)  .<= years[end]) .& .!ismissing(value);

@show  years,sum(sel)

#=
sum(sel)
plot(lon[sel],lat[sel])
plot(lon[sel],lat[sel],"o")



clf();

scatter(lon[sel],lat[sel],10,value[sel]; cmap = "jet"); colorbar()
=#




if false
    fname = "http://sdn.oceanbrowser.net:8081/data/SeaDataNet-domains/Baltic/Salinity.19002012.4Danl.nc"
    varname = "Salinity"
    Slon,Slat,S = surfacemean(fname,varname)

    fname = "http://sdn.oceanbrowser.net:8081/data/SeaDataNet-domains/Baltic/Temperature.19002012.4Danl.nc"
    varname = "Temperature"
    Tlon,Tlat,T = surfacemean(fname,varname)

    ds = Dataset(expanduser("~/workspace/ContourExtractor/data/dist2coast.nc"));
    Dlon = nomissing(ds["lon"][:])
    Dlat = nomissing(ds["lat"][:])

    i0,i1 = extrema(find(Tlon[1]-1 .<= Dlon .<= Tlon[end]+1))
    j0,j1 = extrema(find(Tlat[1]-1 .<= Dlat .<= Tlat[end]+1))

    D = nomissing(ds["distance2coast"][i0:i1,j0:j1]);


    itp = interpolate((Dlon[i0:i1],Dlat[j1:-1:j0]), D[:,end:-1:1], Gridded(Linear()));
    Di = itp[Tlon,Tlat];

    figure(figsize=(11,8))
    subplot(2,2,1)
    pcolor(Tlon,Tlat,T')
    title("Surface mean temperature")
    colorbar()

    subplot(2,2,2)
    pcolor(Slon,Slat,S')
    title("Surface mean salinity")
    colorbar()

    subplot(2,2,3)
    Di[Di .< 0] .= NaN
    pcolor(Tlon,Tlat,Di')
    title("Distance from coast")
    colorbar()


    blon,blat,bath = DIVAnd.load_bath(bathname,bathisglobal,Tlon,Tlat);

    subplot(2,2,4)
    bath[bath .< 0] .= NaN;
    pcolor(Tlon,Tlat,bath')
    title("Depth")
    colorbar()


    MAT.matwrite(joinpath(datadir,"covariables.mat"),
                 Dict("bath" => bath,
                      "T" => T,
                      "S" => S,
                      "Di" => Di,
                      "Tlon" => Tlon,
                      "Tlat" => Tlat))
else
    dd = MAT.matread(joinpath(datadir,"covariables.mat"))

    bath = dd["bath"]
    T = dd["T"]
    S = dd["S"]
    Di = dd["Di"]
    Tlon = dd["Tlon"]
    Tlat = dd["Tlat"]
end

# for plotting
bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,Tlon,Tlat);

dd = matread(joinpath(datadir,"mean2017.mat"))
Clon = dd["lon"];
Clat = dd["lat"];
chl = Float64.(dd["chl"]);

itp = interpolate((Clon,Clat), chl[:,:,1], Gridded(Linear()));
if VERSION >= v"0.7"
    itp = extrapolate(itp,Line());
end
chli = itp(Tlon,Tlat);
minc = minimum(chl[.!isnan.(chl)])
chli[chli .<= minc] .= NaN;


O2 = matread(joinpath(datadir,"covariables_O2_mean.mat"))["O2"]


X,Y = DIVAnd.ndgrid(Tlon,Tlat)

field = cat(T,log.(S),bath,Di,X,Y, dims = 3)
#field = cat(T,S,bath,Di,X,Y, dims = 3)
#field = cat(T,S, dims = 3)

field = cat(T,S,bath,Di,X,Y,log.(chli), dims = 3)
#field = cat(T,S,bath,Di, dims = 3)

#field = cat(T,log.(S),bath,Di,O2,X,Y, dims = 3)


mask = .!isnan.(bath);
mask[1:30,70:end] .= false;
mask[1:4,1:12] .= false;
#mask[:] = true

mask[(19.6 .<= X .<= 20.2) .& (63.4 .<= Y .<= 63.8)] .= true

#pcolor(mask')

fieldf = DIVAnd.ufill(field, .!isnan.(field));

# normalize

for k = 1:size(fieldf,3)
    tmp = fieldf[:,:,k][mask];
    fieldf[:,:,k] = (fieldf[:,:,k] .- mean(tmp))./std(tmp)
end

# 3D extension
mask = repeat(mask,inner=(1,1,length(years)))
fieldf = repeat(reshape(fieldf,size(bath)...,1,size(fieldf,3)),outer=(1,1,length(years),1));

trans,invtrans = DIVAnd.Anam.loglin(30, epsilon = 1);
#trans,invtrans = DIVAnd.Anam.notransform();

if false
    pcolor(fieldf[:,:,1]')



    PyPlot.plt[:hist](value[sel],100)
    title("Histogram");
    savefig(expanduser("~/tmp/Emodnet-Bio/hist_value.png"))



    clf();
    PyPlot.plt[:hist](trans.(value[sel]),100)
    title("Histogram (transformed)");
    savefig(expanduser("~/tmp/Emodnet-Bio/hist_value_trans.png"))

end

mask2,(pm,pn,po),(xi,yi,zi) = DIVAnd.domain(bathname,bathisglobal,Tlon,Tlat,years)

XY = DIVAnd.ndgrid(Tlon,Tlat)

xy = (lon[sel],lat[sel],time[sel]);
v = value[sel];
meanv = binobs(xy,v,XY);



Ntries = 1


RMS_diva = zeros(Ntries)
RMS_diva_covar = zeros(Ntries)



epsilon2 = 0.5

trans_value = trans.(value[sel])
mvalue = mean(trans_value)

varbak,len,dbinfo = DIVAnd.fitlen(
    (lon[sel],lat[sel]),trans_value,0;
    distfun = DIVAnd.distfun_m)

trans_value = trans.(value[sel])
mvalue = mean(trans_value)

f = trans_value .- mvalue

bestfactorl,bestfactore, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D =
    DIVAnd.DIVAnd_cv(mask[:,:,1],(pm[:,:,1],pn[:,:,1]),(xi[:,:,1],yi[:,:,1]),(lon[sel],lat[sel]),f,len,epsilon2,2,3,0);

@warn "fix bestfactore"
bestfactore = 0.2438612845535494

epsilon2 = epsilon2 * bestfactore

if VERSION >= v"0.7"
    Random.seed!(1234)
else
    srand(12345)
end

@show len
@show epsilon2
@show bestfactore

Ntries = 1
value_analysis = zeros(size(mask))
value_analysis2 = zeros(size(mask))

lent = 0.6 # years
niter = 10000
niter = 100000

#for l = 1:Ntries
    l=1

    # 5% for cross-validation
    selcv = rand(Float64,size(sel)) .< 0.10
    #selcv = rand(Float64,size(sel)) .< -0.10
    @show sum(sel .& selcv)

    #    value_analysis,RMS_diva[l],value_analysis2,RMS_diva_covar[l] = comp(selcv)

    lon_a = lon[sel .& .!selcv]
    lat_a = lat[sel .& .!selcv]
    time_a = time[sel .& .!selcv]
    value_a = Vector{Float64}(value[sel .& .!selcv])

    lon_cv = lon[sel .& selcv]
    lat_cv = lat[sel .& selcv]
    time_cv = time[sel .& selcv]
    value_cv = Vector{Float64}(value[sel .& selcv])

    trans_value_a = trans.(value_a)
    mvalue = mean(trans_value_a)

    f = trans_value_a .- mvalue

    value_analysis[:],s = DIVAnd.DIVAndrun(
        mask,(pm,pn,po),(xi,yi,zi),(lon_a,lat_a,time_a),
        f,
        (len,len,lent),epsilon2;
        alphabc = 0,
    );

    value_analysis .= invtrans.(value_analysis .+ mvalue);

    #lon[sel],lat[sel],10,value[sel]
    RMS_diva[l] = validate((Tlon,Tlat,years),value_analysis,
    		           (lon_cv,lat_cv,time_cv),value_cv)

    #epsilon2ap = 0.4
    #epsilon2ap = 1.5
    epsilon2ap = 2
    #epsilon2ap = 5
    #epsilon2ap = 50
    #@show std(value[sel])
    epsilon2ap = epsilon2

    #NLayers = [size(fieldf)[end],3,1]
    NLayers = [size(fieldf)[end],4,1]
    #NLayers = [size(fieldf)[end],5,1]
    #NLayers = [size(fieldf)[end],2,1]
    #NLayers = [size(fieldf)[end],1]
    #NLayers = []
    RMS_diva_covar_his = Float64[]
    stdf = std(f)

    value_analysis2[:],fw0 = analysisprob(
        mask,(pm,pn,po),(xi,yi,zi),(lon_a,lat_a,time_a),
        f / stdf,
        (len,len,lent),epsilon2ap,
        fieldf,
        NLayers,
        costfun = regression,
        niter = niter,
        dropoutprob = 0.,
        L2reg = 0.,
	    plotres = (i,lossi,field,y) -> begin
            #@show i
	        value_analysis2 = invtrans.(stdf * field .+ mvalue)
    	    RMS = validate((Tlon,Tlat,years),value_analysis2,
    	    	           (lon_cv,lat_cv,time_cv),value_cv)
            push!(RMS_diva_covar_his,RMS)
	        @printf("| %10d | %30.5f | %30.5f | %30.5f |\n",i,lossi,RMS,1 - RMS^2/mean(RMS_diva)^2)
	    end,
	    plotevery = 1000
    )

    value_analysis2 .= invtrans.(stdf * value_analysis2 .+ mvalue)

    RMS_diva_covar[l] = validate((Tlon,Tlat,years),value_analysis2,
    		           (lon_cv,lat_cv,time_cv),value_cv)

    #return value_analysis,RMS_diva,value_analysis2,RMS_diva_covar

#end


#value_analysis,RMS_diva_,value_analysis2,RMS_diva_covar_ = comp(falses(size(sel)))


@show mean(RMS_diva), mean(RMS_diva_covar)
@show 1 - mean(RMS_diva_covar)^2/mean(RMS_diva)^2

#=
skill
epsilon2ap = 2
NLayers = [size(fieldf)[end],4,1]
L2reg = 0.,
skill = 0.08603493221507186


epsilon2ap = 1.5
skill = 0.08600723001956245


NLayers = [size(fieldf)[end],3,1]
epsilon2ap = 2
skill = 0.08553042651362863



epsilon2ap = 2
NLayers = [size(fieldf)[end],4,1]
L2reg = 0.,
skill = 0.08592267660472175


with O2
0.08585405392982948
=#

#end
