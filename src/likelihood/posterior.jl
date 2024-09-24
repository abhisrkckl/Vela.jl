export get_lnpost_func

function get_lnpost_func(
    model::TimingModel,
    toas::Vector{T},
    vectorize::Bool = false,
) where {T<:TOABase}
    lnprior_func = get_lnprior_func(model)
    lnlike_func =
        vectorize ? get_lnlike_serial_func(model, toas) : get_lnlike_func(model, toas)
    nfree = model.param_handler._nfree

    function lnpost(params)
        lnpr = lnprior_func(params)
        return isfinite(lnpr) ? lnpr + lnlike_func(params) : -Inf
    end

    if vectorize
        function _lnpost_vector(paramss)
            @assert size(paramss)[2] == nfree
            nparamss = size(paramss)[1]
            result = Vector{Float64}(undef, nparamss)
            @threads for ii = 1:nparamss
                result[ii] = lnpost(paramss[ii, :])
            end
            return result
        end

        return _lnpost_vector
    else
        return lnpost
    end
end
