import DIVAnd

if VERSION >= v"0.7"
    using Test
    using Dates
    using Dates: now
    using Random
else
    using Base.Test
end

using Knet

logistic(x) = 1 / (1 + exp(-x))
logit(x) = log(x / (1 - x))
function weightbias(NLayers,fun = randn)
    # a vector of matrices
    W = Matrix{Float64}[]
    for i = 1:length(NLayers)-1
        push!(W,fun(NLayers[i],NLayers[i+1]))
        push!(W,fun(1,NLayers[i+1]))
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
    J = -sum(y .* log.(HP)  + (1-y) .* log.(1 - HP)) / length(y)
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

        for i = 1:2:length(w)-2
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
    #@show size(H),size(field_test)
    return sum(abs2,H*field_test - y)/epsilon2
end

# likelyhood assuming weight_test is correct
function nll(field_test,H,y,epsilon2)
    prob_test = logistic.(field_test)

    eps = 0.00001
    P_test = (1-2*eps) * prob_test + eps

    #P = DIVAnd.pack(sv,(prob_test,));

    # probability at observed location
    HP = H*P_test

    #@show extrema(HP)
    #@show extrema(prob_test)

    #    δHP = randn(size(HP)) * 1e-6
    #    @show J(y,HP + δHP) - J(y,HP)
    #    @show ∇J(y,HP,δHP)

    return J(y,HP)/epsilon2
end



function background(iB,field_test)
    return field_test' * iB * field_test
end

function ∇background(iB,field_test)
    return 2 * (iB * field_test)
end


function loss_obs(iB,fieldp,H,y,fw,costfun,epsilon2,dropoutprob,L2reg = 0)
    w = fw[1:end-1]
    fieldb = fw[end]

    field_test = model_field(fieldp,fw,dropoutprob)
    #@show size(fieldp),size(field_test)

    J = costfun(field_test,H,y,epsilon2)

    if L2reg != 0
        for i = 1:length(w)
            J += sum(w[i].^2) * L2reg
        end
    end
    return J
end

function loss(iB,fieldp,H,y,fw,costfun,epsilon2,dropoutprob = 0,L2reg = 0)
    J = loss_obs(iB,fieldp,H,y,fw,costfun,epsilon2,dropoutprob,L2reg) + background(iB,fw[end][:,1])
    return J
end


grad_loss_obs = Knet.grad(loss_obs,5)


function grad_loss(iB,fieldp,H,y,fw,costfun,epsilon2,dropoutprob=0,L2reg = 0)
    w = fw[1:end-1]
    fieldb = fw[end]

    δfw = grad_loss_obs(iB,fieldp,H,y,fw,costfun,epsilon2,dropoutprob,L2reg)
    δfw[end] += ∇background(iB,fieldb[:,1])
    return δfw
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

function analysisprob(mask,pmn,xyi,obspos,y,len,epsilon2,field,NLayers;
                      plotres = (i,prob,y) -> nothing,
                      plotevery = -1,
                      niter::Int = 10000,
                      costfun = nll,
                      dropoutprob = 0.,
                      L2reg = 0,
                      learning_rate = 0.01,
                      )

    sv = DIVAnd.statevector((mask,));
    fieldp = DIVAnd.packens(sv,(field,))

    alpha = Float64[]
    moddim = Float64[]
    btrunc = []
    scale_len = true

    # DIVAnd background constraint
    s = DIVAnd.DIVAnd_background(
        Val{:sparse},mask,pmn,len,alpha,moddim,scale_len,[];
        btrunc = btrunc);

    iB = s.iB

    # observation operator
    HI = DIVAnd.localize_separable_grid(obspos,mask,xyi);
    H,out = DIVAnd.sparse_interp(mask,HI);
    H = H * DIVAnd.sparse_pack(mask)'
    # check if out observations are correctly handeld

    @show length(y),sum(out),size(H)

    # random initial weighs and baises with the
    # appropriate size
    weight_bias_test = weightbias(NLayers)
    #weight_bias_test = weightbias(NLayers, (sz...) -> randn(sz...)/100)

    #field_test = DIVAnd.random(mask,pmn,len,1)[:,:,:,1][mask][:,1:1]
    field_test = zeros(sv.n,1)

    fw0 = deepcopy([weight_bias_test..., field_test])

    # @show size(fw0[end])
    # @show loss(iB,fieldp,H,y,fw0,costfun,epsilon2,dropoutprob,L2reg)

    #optim = Knet.optimizers(fw0, Knet.Adam; lr = learning_rate)
    optim = Knet.optimizers(fw0, Knet.Sgd; lr = learning_rate)
    t0 = now()

    fi,s2 = DIVAnd.DIVAndrun(mask,pmn,xyi,obspos,y,len,epsilon2; alphabc = 0)
    # @show sum(s2.obsout)
    x = fw0[end]
    gradloss = grad_loss(iB,fieldp,H,y,fw0,costfun,epsilon2,dropoutprob,L2reg)
    gradloss2 = 2* s2.H' * (s2.R \ (s2.H*x - y)) + 2 * s2.iB * x
    @show gradloss[end][1:10]
    @show gradloss2[1:10]
    @show maximum(abs.(gradloss2 - gradloss[end]))
    # @show gradloss[1:end-1]

    fw0[end] = fi[mask][:,1:1]
    x = fw0[end]
    gradloss = grad_loss(iB,fieldp,H,y,fw0,costfun,epsilon2,dropoutprob,L2reg)
    gradloss2 = 2* s2.H' * (s2.R \ (s2.H*x - y)) + 2 * s2.iB * x
    @show gradloss[end][1:10]
    @show gradloss2[1:10]
    @show maximum(abs.(gradloss2 - gradloss[end]))

    # DIVAnd analysis as a first guess
    fi,s2 = DIVAnd.DIVAndrun(mask,pmn,xyi,obspos,y,len,epsilon2; alphabc = 0)

    @show size.(weight_bias_test)
    weight_bias_test = weightbias(NLayers,zeros);
    @show size.(weight_bias_test)
    @show size.(fw0)

    fw0 = deepcopy([weight_bias_test..., fi[mask][:,1:1]])
    @show size.(fw0)

    for i = 1:niter
        #iobssel = rand(Float64,size(y)) .< 0.1
        iobssel = rand(Float64,size(y)) .< 1

        if ((i-1) % plotevery == 0) && (plotevery != -1)
            prob_estim = model_field(fieldp,fw0,0)
            prob_estim = DIVAnd.unpack(sv,prob_estim)[1]
	        prob_estim[.!mask] .= NaN

            lossi = loss(iB,fieldp,H,y,fw0,costfun,epsilon2,0,L2reg)

            if costfun == nll
                plotres(i-1,lossi,logistic.(prob_estim),y)
            else
                plotres(i-1,lossi,prob_estim,y)
            end
        end


        #gradloss = grad_loss(iB,fieldp,H,y,fw0,costfun,epsilon2,dropoutprob,L2reg)
        gradloss = grad_loss(iB,fieldp,H[iobssel,:],y[iobssel],fw0,costfun,epsilon2,dropoutprob,L2reg)
        @show gradloss[end][1:10]
        @show maximum(gradloss[end])

        x = fw0[end]
        gradloss = [2* s2.H' * (s2.R \ (s2.H*x - y)) + 2 * s2.iB * x]
        @show gradloss[end][1:10]
        @show maximum(gradloss[end])


        lossi = loss(iB,fieldp,H,y,fw0,costfun,epsilon2,0,L2reg)
        @show lossi
        update!(fw0, gradloss, optim)
        lossi = loss(iB,fieldp,H,y,fw0,costfun,epsilon2,0,L2reg)
        @show lossi
        error("ll")
        if (now() - t0) > Dates.Second(3)
            #@show i,loss(iB,fieldp,H,y,fw0,costfun,epsilon2,0,L2reg)
            t0 = now()
        end

    end

    #@show maximum(abs.(fi[mask] - fw0[end]))

    prob_estim = model_field(fieldp,fw0,0)
    prob_estim = DIVAnd.unpack(sv,prob_estim)[1]
    prob_estim[.!mask] .= NaN

    #@show maximum(abs.(fi[mask] - prob_estim[mask]))

    if costfun == nll
        return logistic.(prob_estim),fw0
    else
        return prob_estim,fw0
    end
end
