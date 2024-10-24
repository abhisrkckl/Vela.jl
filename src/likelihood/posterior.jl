export get_lnpost_func


"""
    get_lnpost_func(::TimingModel, toas::Vector{T}, vectorize::Bool = false) where {T<:TOABase}

Returns a callable that evaluates the log-posterior given a collection of parameter values.
If `vectorize` is `true`, then the function supports parallel evaluation on different points
in the parameter space. 
"""
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
            @assert (length(size(paramss)) == 2 && size(paramss)[2] == nfree) "Please pass a 2D array to this function for vectorized evaluation, where each row has `Nfree` elements."
            nparamss = size(paramss)[1]
            result = Vector{Float64}(undef, nparamss)
            @threads :static for ii = 1:nparamss
                result[ii] = lnpost(paramss[ii, :])
            end
            return result
        end

        return _lnpost_vector
    else
        return lnpost
    end
end
