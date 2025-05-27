export Parameter,
    MultiParameter,
    read_param,
    ParamHandler,
    read_params,
    get_free_param_names,
    get_free_param_units,
    get_free_param_labels,
    get_free_param_prefixes,
    read_param_values_to_vector,
    get_scale_factors,
    is_noise,
    get_num_timing_params,
    get_marginalized_param_default_values,
    get_marginalized_param_scale_factors

"""
    Parameter

A single model parameter.

Corresponds to `floatParameter`, `AngleParameter`, or `MJDParameter` in `PINT`.
"""
struct Parameter
    name::Symbol
    default_value::Float64
    dimension::Int
    frozen::Bool
    original_units::String
    unit_conversion_factor::Float64
    noise::Bool
end

Parameter(
    name::Symbol,
    default_quantity::GQ,
    frozen::Bool,
    original_units::String,
    unit_conversion_factor::Float64,
    noise::Bool,
) = Parameter(
    name,
    value(default_quantity),
    udim(default_quantity),
    frozen,
    original_units,
    unit_conversion_factor,
    noise,
)

"""
    MultiParameter

A set of model parameters that are characterized by a common name and a varying index.

Corresponds to `maskParameter` or `prefixParameter` in `PINT`.
"""
struct MultiParameter
    name::Symbol
    parameters::Vector{Parameter}

    function MultiParameter(name::Symbol, parameters::Vector{Parameter})
        @assert allequal([par.noise for par in parameters])
        return new(name, parameters)
    end
end

is_noise(mpar::MultiParameter) = mpar.parameters[1].noise

function _get_single_params_tuple(single_params, noise)
    pkeys = Tuple((sp.name for sp in single_params if sp.noise == noise))
    pquants = Tuple((
        GQ{sp.dimension}(sp.default_value) for sp in single_params if sp.noise == noise
    ))
    return (; zip(pkeys, pquants)...)
end

function _get_multi_params_tuple(multi_params, noise)
    pkeys = Tuple((mp.name for mp in multi_params if is_noise(mp) == noise))
    ptups = Tuple((
        Tuple(_get_single_params_tuple(mp.parameters, noise)) for
        mp in multi_params if is_noise(mp) == noise
    ))
    return (; zip(pkeys, ptups)...)
end

_get_params_tuple(single_params, multi_params) = merge(
    _get_single_params_tuple(single_params, false),
    _get_multi_params_tuple(multi_params, false),
    _get_single_params_tuple(single_params, true),
    _get_multi_params_tuple(multi_params, true),
)

function _get_free_indices(single_params, multi_params)
    idx = 1
    _free_indices = Int[]
    for noise in (false, true)
        @inbounds for spar in single_params
            if spar.noise == noise
                if !spar.frozen
                    push!(_free_indices, idx)
                end
                idx += 1
            end
        end

        @inbounds for mpar in multi_params
            for (jj, param) in enumerate(mpar.parameters)
                if param.noise == noise
                    if !param.frozen
                        push!(_free_indices, idx)
                    end
                    idx += 1
                end
            end
        end
    end
    return _free_indices
end

"""
    ParamHandler{PT<:NamedTuple}

Handles the creation of a parameter tuple from a collection of free parameter values.
"""
struct ParamHandler{ParamsType<:NamedTuple}
    single_params::Vector{Parameter}
    multi_params::Vector{MultiParameter}
    _default_params_tuple::ParamsType
    _default_values::Vector{Float64}
    _free_indices::Vector{Int}
    _nfree::Int
end

function ParamHandler(single_params, multi_params)
    default_params = _get_params_tuple(single_params, multi_params)
    default_values = collect(map(value, Iterators.flatten(default_params)))

    free_indices = _get_free_indices(single_params, multi_params)

    ParamHandler(
        single_params,
        multi_params,
        default_params,
        default_values,
        free_indices,
        length(free_indices),
    )
end

"""
    read_params(ph::ParamHandler{PT}, free_values)::PT where {PT<:NamedTuple}

Create a parameter tuple from a collection of free parameter values.

Reverse of `read_param_values_to_vector()`.
"""
function read_params(
    ph::ParamHandler{ParamsType},
    free_values,
)::ParamsType where {ParamsType<:NamedTuple}
    @assert length(free_values) == ph._nfree
    values = copy(ph._default_values)
    @inbounds for (idx, value) in zip(ph._free_indices, free_values)
        values[idx] = value
    end

    return reinterpret(ParamsType, values)[1]
end

function _get_free_param_attribute(param_handler::ParamHandler, attr_func)
    pattrs = []

    for noise in (false, true)
        @inbounds for spar in param_handler.single_params
            if !spar.frozen && spar.noise == noise
                push!(pattrs, attr_func(spar, spar))
            end
        end

        @inbounds for mpar in param_handler.multi_params
            for param in mpar.parameters
                if !param.frozen && param.noise == noise
                    push!(pattrs, attr_func(mpar, param))
                end
            end
        end
    end

    return pattrs
end

"""
    get_free_param_names(param_handler::ParamHandler)::Vector{String}

Generate an ordered collection of free parameter names."""
get_free_param_names(param_handler::ParamHandler)::Vector{String} =
    _get_free_param_attribute(param_handler, (mpar, param) -> string(param.name))

"""
    get_free_param_prefixes(param_handler::ParamHandler)::Vector{String}

Generate an ordered collection of free parameter prefixes."""
get_free_param_prefixes(param_handler::ParamHandler)::Vector{String} =
    _get_free_param_attribute(param_handler, (mpar, param) -> string(mpar.name))

"""
    get_free_param_labels(param_handler::ParamHandler; delim::String = "\n")::Vector{String}

Generate an ordered collection of free parameter labels."""
get_free_param_labels(param_handler::ParamHandler; delim::String = "\n")::Vector{String} =
    _get_free_param_attribute(
        param_handler,
        (mpar, param) -> (
            isempty(param.original_units) || param.original_units=="1" ?
            string(param.name) : "$(string(param.name))$delim($(param.original_units))"
        ),
    )

"""
    get_free_param_units(param_handler::ParamHandler)::Vector{String}

Generate an ordered collection of free parameter units (astropy-compatible)."""
get_free_param_units(param_handler::ParamHandler)::Vector{String} =
    _get_free_param_attribute(param_handler, (mpar, param) -> param.original_units)

"""
    read_param_values_to_vector(param_handler::ParamHandler, params::NamedTuple)::Vector{Float64}

Generate a collection of free parameter values from a parameter tuple.

Reverse of `read_params()`
"""
function read_param_values_to_vector(
    param_handler::ParamHandler,
    params::NamedTuple,
)::Vector{Float64}
    param_vec = Float64[]

    for noise in (false, true)
        @inbounds for spar in param_handler.single_params
            if !spar.frozen && spar.noise == noise
                push!(param_vec, Float64(params[spar.name].x))
            end
        end

        @inbounds for mpar in param_handler.multi_params
            for (jj, param) in enumerate(mpar.parameters)
                if !param.frozen && param.noise == noise
                    push!(param_vec, Float64(params[mpar.name][jj].x))
                end
            end
        end
    end

    return param_vec
end
read_param_values_to_vector(param_handler::ParamHandler) =
    read_param_values_to_vector(param_handler, param_handler._default_params_tuple)

"""
    get_scale_factors(param_handler::ParamHandler)::Vector{Float64}

Get the scale factors that convert the free parameters from `Vela.jl`'s 
internal representation to the units used in `PINT`."""
get_scale_factors(param_handler::ParamHandler)::Vector{Float64} =
    _get_free_param_attribute(param_handler, (mpar, param) -> param.unit_conversion_factor)

function get_num_timing_params(param_handler::ParamHandler)
    ntmdim = 0

    @inbounds for spar in param_handler.single_params
        if !spar.frozen && !spar.noise
            ntmdim += 1
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen && !param.noise
                ntmdim += 1
            end
        end
    end

    return ntmdim
end

function get_marginalized_param_default_values(
    param_handler::ParamHandler,
    marginalized_param_names::Vector{String},
)
    result = zeros(Float64, length(marginalized_param_names))

    @inbounds for spar in param_handler.single_params
        pname = string(spar.name)
        idx = findfirst(x -> x==pname, marginalized_param_names)
        if !isnothing(idx)
            result[idx] = spar.default_value
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            pname = string(param.name)
            idx = findfirst(x -> x==pname, marginalized_param_names)
            if !isnothing(idx)
                result[idx] = param.default_value
            end
        end
    end

    return result
end

function get_marginalized_param_scale_factors(
    param_handler::ParamHandler,
    marginalized_param_names::Vector{String},
)
    result = zeros(Float64, length(marginalized_param_names))

    @inbounds for spar in param_handler.single_params
        pname = string(spar.name)
        idx = findfirst(x -> x==pname, marginalized_param_names)
        if !isnothing(idx)
            result[idx] = spar.unit_conversion_factor
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            pname = string(param.name)
            idx = findfirst(x -> x==pname, marginalized_param_names)
            if !isnothing(idx)
                result[idx] = param.unit_conversion_factor
            end
        end
    end

    return result
end
