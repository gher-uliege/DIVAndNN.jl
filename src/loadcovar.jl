function loadcovar(xyi,covars_fname;
                   covars_coord = false,
                   covars_const = false)

    gridlon,gridlat = xyi[1:2]
    ndimensions = length(xyi)
    domain_size = length.(xyi)

    ncovars = length(covars_fname)

    if covars_const
        ncovars = ncovars + 1
    end

    if covars_coord
        ncovars = ncovars + ndimensions
    end

    sz = (domain_size...,ncovars)

    field = zeros(sz)

    for i = 1:length(covars_fname)
        (fname,varname,trans) = covars_fname[i]

        Dataset(fname) do ds
            tmp = nomissing(ds[varname][:],NaN)
            tmp = trans.(tmp)
            if ndimensions == 3
                if ndims(tmp) == 2
                    field[:,:,:,i] = repeat(tmp,inner = (1,1,domain_size[3]))
                else
                    field[:,:,:,i] = tmp
                end
            else
                field[:,:,i] = tmp
            end
        end
    end

    if covars_const
        i = length(covars_fname)+1
        @info "add const covariable"

        if ndimensions == 3
            field[:,:,:,i] .= 1
        else
            field[:,:,i] .= 1
        end
    end

    if covars_coord
        XYI = DIVAnd.ndgrid(xyi...)

        if ndimensions == 3
            @info "add lon/lat/time as covariable"
            field[:,:,:,end-2] = XYI[1]
            field[:,:,:,end-1] = XYI[2]
            field[:,:,:,end]   = XYI[3]
        else
            @info "add lon/lat as covariable"
            field[:,:,end-1] = XYI[1]
            field[:,:,end] = XYI[2]
        end
    end

    return field
end
