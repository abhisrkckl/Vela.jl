struct NonCentralParametrization{GPComponentsTuple<:Tuple,IdxInfo<:NamedTuple}
    gp_components::GPComponentsTuple
    M::Matrix{Float64}
    Phidiag::Vector{Float64}
    y::Vector{Float64}
    a0::Vector{Float64}
    gp_amplitude_index_info::IdxInfo
    ntmpar::Int
end

function NonCentralParametrization(
    model::TimingModel,
    M::Matrix,
    y::Vector,
    timing_param_weight::Float64,
)
    gp_components = filter(is_gp_noise, model.components)

    param_prefixes = get_free_param_prefixes(model)

    index_info_keys = Symbol[]
    index_info_vals = NTuple{2,Int}[]
    for gpcomp in gp_components
        ampl_params = gp_ampl_params(gpcomp)
        for ampl_par in ampl_params
            ampl_par_str = string(ampl_par)
            start_idx = findfirst(==(ampl_par_str), param_prefixes)
            nparams = count(==(ampl_par_str), param_prefixes)
            append!(index_info_keys, ampl_par)
            append!(index_info_vals, (start_idx, nparams))
        end
    end
    index_info = (; zip(index_info_keys, index_info_vals)...)

    ntmpar = get_num_timing_params(model)
    Phidiag = fill(timing_param_weight, ntmpar)
    for (start_idx, npar) in index_info_vals
        Phidiag[start_idx:(start_idx+npar)] = 1.0
    end

    a0 = read_param_values_to_vector(model)[1:ntmpar]

    @assert size(M)[1] == length(y) "The design matrix shape does not match the number of residuals."
    @assert size(M)[2] == ntmpar "The design matrix shape does not match the number of timing parameters."

    return NonCentralParametrization(gp_components, M, Phidiag, y, a0, index_info, ntmpar)
end

function _scale_designmatrix!(M::Matrix, weights::Vector, start_idx, nparams)
    Ntim, Npar = size(M)
    @assert length(weights) == nparams
    @assert start_idx >= 1 && start_idx + nparams <= Npar

    @inbounds for (p, weight) in enumerate(weights)
        @simd for j = 1:Ntim
            M[j, start_idx+p] *= weight
        end
    end
end

function designmatrix(ncp::NonCentralParametrization, params::NamedTuple)::Matrix
    M = copy(ncp.M)

    for gpcomp in ncp.gp_components
        ampl_params = gp_ampl_params(gpcomp)
        weightss = gp_weights(gpcomp, params)

        for (par, weights) in zip(ampl_params, weightss)
            start_idx, nparams = ncp.gp_amplitude_index_info[par]
            _scale_designmatrix!(M, weights, start_idx, nparams)
        end
    end

    return M
end

function apply_noncentral_transform(
    model::TimingModel,
    ncp::NonCentralParametrization,
    params,
)
    alpha = params[1:ncp.ntmpar]

    params_0 = read_params(model, params) # Only noise parameters are valid here
    Phidiag = ncp.Phidiag
    M = designmatrix(ncp, params_0)
    y = ncp.y
    a0 = ncp.a0

    return a0 + calc_noncentral_transform(M, Ndiag, Phidiag, y, alpha)
end

function get_noncentral_param_initial_sample(model::TimingModel)
    ntmpar = get_num_timing_params(model)
    cube = rand(model.param_handler._nfree)
    sample = prior_transform(model, cube)
    sample[:ntmpar] = randn(ntmpar)
    return sample
end
