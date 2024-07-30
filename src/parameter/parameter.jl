export Parameter,
    MultiParameter,
    read_param,
    ParamHandler,
    read_params,
    get_free_param_names,
    read_param_values_to_vector,
    get_scale_factors

"""A single model parameter.

Corresponds to `floatParameter`, `AngleParameter`, or `MJDParameter` in `PINT`.
"""
struct Parameter
    name::Symbol
    default_quantity::GQ{Float64}
    frozen::Bool
    original_units::String
    unit_conversion_factor::Float64
end

"""A set of model parameters that are characterized by a common name and a varying index.

Corresponds to `maskParameter` or `prefixParameter` in `PINT`.
"""
struct MultiParameter
    name::Symbol
    parameters::Vector{Parameter}
end

"""Read a parameter value into a `GQ` object."""
function read_param(param::Parameter, value::Float64)::GQ{Float64}
    @assert !param.frozen "Refusing to read a frozen parameter $(param.name)."

    return quantity_like(param.default_quantity, value)
end

function _get_single_params_tuple(single_params)
    pkeys = Tuple((sp.name for sp in single_params))
    pquants = Tuple((sp.default_quantity for sp in single_params))
    return (; zip(pkeys, pquants)...)
end

function _get_multi_params_tuple(multi_params)
    pkeys = Tuple((mp.name for mp in multi_params))
    ptups = Tuple((Tuple(_get_single_params_tuple(mp.parameters)) for mp in multi_params))
    return (; zip(pkeys, ptups)...)
end

_get_params_tuple(single_params, multi_params) =
    merge(_get_single_params_tuple(single_params), _get_multi_params_tuple(multi_params))

"""Handles the creation of a parameter tuple from a collection of free parameter values."""
struct ParamHandler{ParamsType<:NamedTuple}
    single_params::Vector{Parameter}
    multi_params::Vector{MultiParameter}
    _default_params_tuple::ParamsType
    _default_quantities::Vector{GQ{Float64}}
    _free_indices::Vector{Int}
    _nfree::Int
end

function ParamHandler(single_params, multi_params)
    default_params = _get_params_tuple(single_params, multi_params)
    default_quantities = collect(Iterators.flatten(default_params))

    all_params = [
        single_params
        collect(Iterators.flatten([mpar.parameters for mpar in multi_params]))
    ]
    free_indices = collect(1:length(all_params))[[!par.frozen for par in all_params]]

    ParamHandler(
        single_params,
        multi_params,
        default_params,
        default_quantities,
        free_indices,
        length(free_indices),
    )
end

"""Create a parameter tuple from a collection of free parameter values.

Reverse of `read_param_values_to_vector`.
"""
function read_params(
    ph::ParamHandler{ParamsType},
    free_values,
)::ParamsType where {ParamsType<:NamedTuple}
    @assert length(free_values) == ph._nfree
    quantities = copy(ph._default_quantities)
    for (idx, value) in zip(ph._free_indices, free_values)
        d = quantities[idx].d
        quantities[idx] = GQ(value, d)
    end

    return reinterpret(ParamsType, quantities)[1]
end

"""Generate an ordered collection of free parameter names."""
function get_free_param_names(param_handler::ParamHandler)::Vector{String}
    pnames = Vector{String}()

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(pnames, string(spar.name))
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                push!(pnames, string(param.name))
            end
        end
    end

    return pnames
end

"""Generate a collection of free parameter values from a parameter tuple.

Reverse of `read_params`
"""
function read_param_values_to_vector(
    param_handler::ParamHandler,
    params::NamedTuple,
)::Vector{Float64}
    param_vec = Float64[]

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(param_vec, Float64(params[spar.name].x))
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for (jj, param) in enumerate(mpar.parameters)
            if !param.frozen
                push!(param_vec, Float64(params[mpar.name][jj].x))
            end
        end
    end

    return param_vec
end
read_param_values_to_vector(param_handler::ParamHandler) =
    read_param_values_to_vector(param_handler, param_handler._default_params_tuple)

function get_scale_factors(param_handler::ParamHandler)
    scale_factors = Float64[]

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(scale_factors, spar.unit_conversion_factor)
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                push!(scale_factors, param.unit_conversion_factor)
            end
        end
    end

    return scale_factors
end
