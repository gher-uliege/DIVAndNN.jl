using DIVAnd
using Test
using DIVAndNN


# grid of background field
mask, pmn, xyi = DIVAnd_squaredom(2, range(-1, stop = 1, length = 50))

n = 100
xy = (zeros(n),zeros(n))
f = zeros(n)
f[1:(n÷4)] .= 1

# correlation length
len = (0.15,0.15)
epsilon2 = 1.

# no covariables
field = zeros(size(mask)...,0)

NLayers = []

value_analysis,fw0 = DIVAndNN.analysisprob(
    mask,pmn,xyi,xy,
    f,
    len,epsilon2,
    field,
    NLayers,
    costfun = DIVAndNN.nll,
    niter = 1000,
    learning_rate = 0.001,
    rmaverage = true,
)

@test minimum(value_analysis) ≈ 1/4 atol=2e-3
@test maximum(value_analysis) ≈ 1/4 atol=2e-3



