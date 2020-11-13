import DIVAnd

using Test
using Dates
using Dates: now
using Random
using Statistics
using Knet

logistic(x) = 1 / (1 + exp(-x))
logit(x) = log(x / (1 - x))
function weightbias(NLayers,fun = randn)
    # a vector of matrices
    W = Matrix{Float64}[]
    for i = 1:length(NLayers)-1
        push!(W,fun((NLayers[i],NLayers[i+1])))
        push!(W,fun((1,NLayers[i+1])))
    end
    return W
end

# negative log likelyhood
# -log of a Binomial distribution

# can be optimized since y_i is just 0 and 1
#=
J = -\frac{1}{m} \sum_i y_i log(p_i) + (1 - y_i) log(1 - p_i)
=#
function J(y,HP)
    J = -sum(y .* log.(HP)  + (1 .- y) .* log.(1 .- HP)) / length(y)
    return J
end

#=
δJ = - \frac{1}{m} \sum_i (\frac{ y_i}{ p_i}  - \frac{1 - y_i}{1 - p_i} ) δp
=#

function ∇J(y,HP,δHP)
    J = -sum( (y ./ HP  - (1-y) ./ (1 - HP) ) .* δHP) / length(y)
    return J
end

function adj_J(y,HP,∂J)
    δHP = - (y ./ HP  - (1-y) ./ (1 - HP) * ∂J ) / length(y)
    return δHP
end

#const epsback = 0.001

relu(x) = max(0.,x)

"""
    fieldp: packed fields of covariables (Nxy x Nfield)
    fw[1:end-1]: neural network weight and biases
    fw[end]: background field (Nxy x 1)
"""
function model_field(fieldp,fw, dropoutprob = 0.)
    w = fw[1:end-1]
    fieldb = fw[end]

    if length(w) > 0
        field_test = fieldp

        for i = 1:2:(length(w)-2)
            field_test = relu.((field_test * w[i] .+ w[i+1]))
            Knet.dropout(field_test, dropoutprob)
        end
        field_test = (field_test * w[end-1] .+ w[end])

        #@show size(DIVAnd.unpackens(s.sv,field_test)[1])
        field_test += fieldb
        return field_test[:,1]
    else
        return fieldb[:,1]
    end
end


function regression(field_test,H,y,epsilon2)
    J = sum(abs2,H*field_test - y)/epsilon2
    #@show size(H),size(field_test)
    return J
end

# likelyhood assuming weight_test is correct
function nll(field_test,H,y,epsilon2)
    #@show extrema(field_test)
    prob_test = logistic.(field_test)

    #eps = 0.00001
    eps = 0.0001
    P_test = (1-2*eps) * prob_test .+ eps
    #P_test = max.(min.(prob_test,1 - eps),eps)

    #P = DIVAnd.pack(sv,(prob_test,));

    # probability at observed location
    HP = H*P_test

    #@show extrema(P_test)
    #@show extrema(HP)
    #@show extrema(prob_test)

    #    δHP = randn(size(HP)) * 1e-6
    #    @show J(y,HP + δHP) - J(y,HP)
    #    @show ∇J(y,HP,δHP)

    return J(y,HP)/epsilon2
end



function background(iB,field_test,field_background)
    return (field_test - field_background)' * iB * (field_test - field_background)
end

function ∇background(iB,field_test,field_background)
    return 2 * (iB * (field_test - field_background))
end


function loss_obs(iB,fieldp,H,y,fw,costfun,epsilon2,dropoutprob,L2reg = 0)
    w = fw[1:end-1]
    fieldb = fw[end]

    #@show dropoutprob
    field_test = model_field(fieldp,fw,dropoutprob)
    ##@show size(fieldp),size(field_test)
    #@show extrema(field_test)

    J = costfun(field_test,H,y,epsilon2)

    if L2reg != 0
        for i = 1:length(w)
            J += sum(w[i].^2) * L2reg
        end
    end
    return J
end

function loss(iB,fieldp,H,y,fw,costfun,epsilon2,field_background,dropoutprob = 0,L2reg = 0)
    Jobs = loss_obs(iB,fieldp,H,y,fw,costfun,epsilon2,dropoutprob,L2reg)
    Jbackground = background(iB,fw[end][:,1],field_background)

    J = Jobs + Jbackground
    return J
end


grad_loss_obs = Knet.grad(loss_obs,5)


function grad_loss(iB,fieldp,H,y,fw,costfun,epsilon2,field_background,dropoutprob=0,L2reg = 0)
    w = fw[1:end-1]
    fieldb = fw[end]

    δfw = grad_loss_obs(iB,fieldp,H,y,fw,costfun,epsilon2,dropoutprob,L2reg)

    δfw[end] += ∇background(iB,fieldb[:,1],field_background)
    return δfw

    #δfw2 = vcat(δfw[1:end-1],[δfw[end] + ∇background(iB,fieldb[:,1])])
    #return δfw2
end

function sampleobs(Nobs,xyi,prob)
    n = ndims(prob)
    sz = size(prob)
    obspos = ntuple(i-> zeros(Nobs),n)
    y = zeros(Int,Nobs)

    index = zeros(Int,n)

    for i = 1:Nobs
        index .= Int[rand(1:sz[i]) for i in 1:n];

        for j = 1:n
            obspos[j][i] = xyi[j][index...]
        end

        y[i] = rand() < prob[index...]
    end

    return obspos,y
end


function sampleobsfield(Nobs,xyi,field)
    n = ndims(xyi[1])
    sz = size(xyi[1])
    obspos = ntuple(i-> zeros(Nobs),n)
    y = zeros(Float64,Nobs)

    index = zeros(Int,n)

    for i = 1:Nobs
        index .= Int[rand(1:sz[i]) for i in 1:n];

        for j = 1:n
            obspos[j][i] = xyi[j][index...]
        end

        y[i] = field[index...]
    end

    return obspos,y
end


"""
    prob_estim,fw0 = DIVAndNN.analysisprob(mask,pmn,xyi,obspos,y,len,epsilon2,field,NLayers;...)

The parameters `mask`,`pmn`,`xyi`,`obspos`,`y`,`len`,`epsilon2` are the same than for
DIVAnd.DIVAndrun. `field` is an array of covariable with the same size of mask but with an extra dimension.
For example if mask has the size 100x101 and if there are 3 co-variables, then field would have the size 100x101x3.
`NLayers` is the number of layer in the neural network. `NLayers[1]` must be equal tot he number of co-variables
and `NLayers[end]` must be 1.

## Keyword parameters

* `niter` (default 10000): number of iterations
* `costfun` (default `nll`): `nll` for probabilities (between 0 and 1) and `regression` for unbounded variables
* `dropoutprob` (default 0.): drop-out probability to avoid overfitting
* `L2reg` (default 0): L2 regularization
* `learning_rate` (default 0.01)
* `maxgrad` (default 5000.): maximum gradient (in absolute values) to stabilize the learning
* `rmaverage` (default false): remove average
* `trainfrac` (default 0.1): fraction of training data
* `epsilon2_background` (default 0.1): `epsilon2` parameter for initial analysis
* `plotres`: is a function that is executed every `plotevery` times to minitor the convergence with the parameters:
   `i`,`lossi`,`prob_estim`,`y`,`gradloss`,`out`,`iobssel`,`obspos`

"""
function analysisprob(mask,pmn,xyi,obspos,y,len,epsilon2,field,NLayers;
                      plotres = (i,prob,y,params...) -> nothing,
                      plotevery = -1,
                      niter::Int = 10000,
                      costfun = nll,
                      dropoutprob = 0.,
                      L2reg = 0,
                      learning_rate = 0.01,
                      maxgrad = 5000.,
                      rmaverage = false,
                      trainfrac = 0.1,
                      epsilon2_background = 0.1,
                      )

    sv = DIVAnd.statevector((mask,));
    fieldp = DIVAnd.packens(sv,(field,))

    alpha = Float64[]
    alpha = Float64[0,2,1]
    moddim = Float64[]
    btrunc = []
    scale_len = true

    # DIVAnd background constraint
    s = DIVAnd.DIVAnd_background(
        Val{:sparse},mask,pmn,len,alpha,moddim,scale_len,[];
        btrunc = btrunc);

    iB = s.iB

    #@show extrema(iB*ones(size(iB,1)))

    # observation operator
    HI = DIVAnd.localize_separable_grid(obspos,mask,xyi);
    H,out = DIVAnd.sparse_interp(mask,HI);
    H = H * DIVAnd.sparse_pack(mask)'
    # check if out observations are correctly handeld

    #@show length(y),sum(out),size(H)
    # remove obs out of grid
    y = y[.!out]
    H = H[.!out,:]
    obspos = map(p -> p[.!out],obspos)



    # DIVAnd analysis as a first guess
    meany = mean(y)
    if rmaverage
        field_background = fill(meany,sv.n)
        ya = y .- meany
    else
        field_background = fill(0.,sv.n)
        ya = y
    end

    fi,s2 = DIVAnd.DIVAndrun(mask,pmn,xyi,obspos,ya,len,epsilon2_background;
                             alpha = alpha,
                             alphabc = 0)
    if rmaverage
        fi = fi .+ meany
    end

    if costfun == nll
        eps = 0.0001
        fi = logit.(clamp.(fi,eps,1-eps))
    end

    # random initial weighs and baises with the
    # appropriate size
    weight_bias_test = weightbias(NLayers,sz -> 0.0001*randn(sz));
    fw0 = deepcopy([weight_bias_test..., fi[mask][:,1:1]])
    optim = Knet.optimizers(fw0, Knet.Adam; lr = learning_rate)

    #@show size.(fw0)

    for i = 1:niter
        iobssel = rand(Float64,size(y)) .<= trainfrac

        gradloss = grad_loss(iB,fieldp,H[iobssel,:],y[iobssel],fw0,costfun,
                             epsilon2,field_background,dropoutprob,L2reg)

        if ((i-1) % plotevery == 0) && (plotevery != -1)
            prob_estim = model_field(fieldp,fw0,0)
            prob_estim = DIVAnd.unpack(sv,prob_estim)[1]
	        prob_estim[.!mask] .= NaN

            #@show extrema(fieldp),extrema(y)
            #@show extrema.(fw0)
            lossi = loss(iB,fieldp,H,y,fw0,costfun,epsilon2,field_background,0,L2reg)
            #@show i,lossi
            if costfun == nll
                plotres(i-1,lossi,logistic.(prob_estim),y,gradloss,out,iobssel,obspos)
            else
                plotres(i-1,lossi,prob_estim,y,gradloss,out,iobssel,obspos)
            end
        end



        for (i,gl) in enumerate(gradloss)
            #@show i,extrema(gl),size(gl)
            clamp!(gl,-maxgrad,maxgrad)
            #@show i,extrema(gl),size(gl)
        end

        update!(fw0, gradloss, optim)
    end

    #@show maximum(abs.(fi[mask] - fw0[end]))

    prob_estim = model_field(fieldp,fw0,0)
    prob_estim = DIVAnd.unpack(sv,prob_estim)[1]
    prob_estim[.!mask] .= NaN

    #@show sum(prob_estim)
    #@show maximum(abs.(fi[mask] - prob_estim[mask]))

    if costfun == nll
        return logistic.(prob_estim),fw0
    else
        return prob_estim,fw0
    end
end
