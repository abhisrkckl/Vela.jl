export calc_lnpost,
    calc_lnpost_serial,
    calc_lnpost_vectorized,
    get_lnpost_func,
    calc_lnpost_transformed,
    calc_lnpost_transformed_vectorized


function calc_lnpost(model::TimingModel, toas::Vector{T}, params) where {T<:TOABase}
    lnpr = calc_lnprior(model, params)
    return !isfinite(lnpr) ? lnpr : lnpr + calc_lnlike(model, toas, params)
end

function calc_lnpost_serial(model::TimingModel, toas::Vector{T}, params) where {T<:TOABase}
    lnpr = calc_lnprior(model, params)
    return !isfinite(lnpr) ? lnpr : lnpr + calc_lnlike_serial(model, toas, params)
end

function calc_lnpost_vectorized(
    model::TimingModel,
    toas::Vector{T},
    paramss,
) where {T<:TOABase}
    nparamss = size(paramss)[1]
    result = Vector{Float64}(undef, nparamss)
    @threads :static for ii = 1:nparamss
        result[ii] = calc_lnpost_serial(model, toas, paramss[ii, :])
    end # COV_EXCL_LINE
    return result
end

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
    return if vectorize
        paramss -> calc_lnpost_vectorized(model, toas, paramss)
    elseif nthreads() == 1
        params -> calc_lnpost_serial(model, toas, params)
    else
        params -> calc_lnpost(model, toas, params) # COV_EXCL_LINE
    end
end

"""
    calc_lnpost_transformed(model::TimingModel, toas::Vector{T}, cube) where {T<:TOABase}

Compute the prior-transformed log-posterior on a sample taken from the unit hypercube. 
"""
function calc_lnpost_transformed(
    model::TimingModel,
    toas::Vector{T},
    cube,
) where {T<:TOABase}
    if !all(x -> 0<=x<=1, cube)
        return -Inf
    end
    params = prior_transform(model, cube)
    return calc_lnlike_serial(model, toas, params)
end

function calc_lnpost_transformed_vectorized(
    model::TimingModel,
    toas::Vector{T},
    cubes,
) where {T<:TOABase}
    ncubes = size(cubes)[1]
    result = Vector{Float64}(undef, ncubes)
    @threads :static for ii = 1:ncubes
        result[ii] = calc_lnpost_transformed(model, toas, cubes[ii, :])
    end # COV_EXCL_LINE
    return result
end
