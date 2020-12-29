function std_or_1(tmp)
    s = std(tmp)
    if s == 0
        return one(s)
    else
        return s
    end
end

"""
    normalize!(field)

"""
function normalize!(mask,field)
    for n = 1:size(field)[end]
        if ndims(field) == 4
            tmp = field[:,:,:,n][mask]
            field[:,:,:,n] = (field[:,:,:,n] .- mean(tmp)) ./ std_or_1(tmp)
        else
            tmp = field[:,:,n][mask];
            field[:,:,n] = (field[:,:,n] .- mean(tmp)) ./ std_or_1(tmp)
        end
    end
end
