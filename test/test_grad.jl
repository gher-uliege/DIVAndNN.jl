using DIVAnd
using Test
using DIVAndNN
using Statistics

# grid of background field
mask, pmn, xyi = DIVAnd_squaredom(2, range(-1, stop = 1, length = 50))

n = 100
xy = (rand(n),rand(n))
f = rand(n)

# correlation length
len = (0.3,0.3)
epsilon2 = 1.

# no covariables
field = zeros(size(mask)...,0)

NLayers = []

gradloss_iter = []
function plotres2(i,lossi,value_analysis,y,gradloss,out,iobssel,obspos)
    push!(gradloss_iter,gradloss)
end

value_analysis,fw0 = DIVAndNN.analysisprob(
    mask,pmn,xyi,xy,
    f,
    len,epsilon2,
    field,
    NLayers,
    costfun = DIVAndNN.regression,
    niter = 1,
    learning_rate = 0.001,
    rmaverage = true,
    epsilon2_background = epsilon2,
    plotres = plotres2,
    plotevery = 1,
    trainfrac = 1,
)

@test maximum(abs.(gradloss_iter[1][end])) ≈ 0 atol=1e-8

obspos = xy
y = f
costfun = DIVAndNN.regression
dropoutprob = 0.
L2reg = 0.

sv = DIVAnd.statevector((mask,));
fieldp = DIVAnd.packens(sv,(field,))

x = DIVAnd.random(mask,pmn,len,1)[:,:,1][mask][:,1:1]
fw0 = [x]

#figure();pcolor(DIVAnd.unpack(sv,x[:,1])[1])

# DIVAnd analysis as a first guess
meany = mean(y)
field_background = fill(meany,sv.n)
ya = y .- meany

fi,s2 = DIVAnd.DIVAndrun(mask,pmn,xyi,obspos,ya,len,epsilon2;
                         alphabc = 0)

gradloss = DIVAndNN.grad_loss(s2.iB,fieldp,s2.H,y,fw0,costfun,epsilon2,field_background,dropoutprob,L2reg)
gradloss2 = 2 * s2.H' * (s2.R \ (s2.H*x - y)) + 2 * s2.iB * (x - field_background)

@test maximum(abs.(gradloss2 - gradloss[end])) ≈ 0 atol=1e-10

